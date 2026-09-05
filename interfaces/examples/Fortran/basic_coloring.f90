! basic_coloring.f90 — minimal example: color the columns of a sparse matrix.
!
! A is the 6x6 tridiagonal matrix tridiag(-1, 2, -1).  Two columns may share a
! color only if they have no nonzero in a common row, so that the corresponding
! columns of a Jacobian can be recovered from a single directional derivative.
! Only the *sparsity pattern* is needed here: smc_coloring never looks at the
! numerical values (see compress_decompress.f90 for the values).
!
! This is the Fortran counterpart of examples/C/basic_coloring.c, with one
! deliberate difference: it sets opts%index_base = 1, so that colptr, rowval and
! the group members returned by the queries are all 1-based — Fortran's own
! convention.  Color *labels* are never shifted by index_base; only indices are.
!
! Compile (after building libsmc with juliac — see interfaces/README.md):
!
!   gfortran -O2 -o interfaces/build/basic_coloring_f \
!       interfaces/examples/Fortran/basic_coloring.f90 \
!       -I interfaces/include \
!       interfaces/build/lib/libsmc.so \
!       -Wl,-rpath,'$ORIGIN/lib'
!
! -I points gfortran at the directory holding smc.f90, the include file below.
! The rpath makes the executable find the bundle without LD_LIBRARY_PATH: it
! assumes the binary sits in interfaces/build/, next to lib/.  The library then
! locates the embedded Julia runtime through its own runpath.  On macOS use
! -Wl,-rpath,@loader_path/lib instead ($ORIGIN is a Linux spelling).
!
! Expected output:
!   SparseMatrixColorings 0.4.x
!   ncolors = 3
!   column colors = [ 1 2 3 1 2 3 ]
!   group 1 = { 1 4 }
!   group 2 = { 2 5 }
!   group 3 = { 3 6 }
!
! Exit code: 0 on success, 1 if any call fails.

program basic_coloring
  use iso_c_binding
  use iso_fortran_env, only: error_unit
  implicit none
  include 'smc.f90'          ! <- after implicit none; brings in the interfaces,
                             !    the enumerators and type(SmcColoringOptions)

  ! -------------------------------------------------------------------------
  ! Problem data: A = tridiag(-1, 2, -1), 6x6, in compressed sparse column form.
  !
  ! With opts%index_base == 1 both colptr and rowval are 1-based, exactly like a
  ! Fortran array:
  !
  !   column j occupies the entries colptr(j) .. colptr(j+1) - 1
  !   rowval(k) is the 1-based row index of entry k
  !
  ! Column 1 touches rows {1,2}, column 2 rows {1,2,3}, and so on.  (The C
  ! example uses index_base = 0, the default, and subtracts one everywhere.)
  !
  ! Every array whose address is handed to C through c_loc must carry the
  ! target attribute — that is what makes taking its address legal, and what
  ! stops the compiler from keeping it in a register or copying it.
  ! -------------------------------------------------------------------------

  integer(c_int), parameter :: m   = 6      ! number of rows
  integer(c_int), parameter :: n   = 6      ! number of columns
  integer(c_int), parameter :: nnz = 16     ! number of stored entries

  integer(c_int), target :: colptr(n+1) = [ 1, 3, 6, 9, 12, 15, 17 ]
  integer(c_int), target :: rowval(nnz) = [ 1, 2,          &
                                            1, 2, 3,       &
                                            2, 3, 4,       &
                                            3, 4, 5,       &
                                            4, 5, 6,       &
                                            5, 6 ]

  type(SmcColoringOptions), target :: opts

  ! Scalars written by the library through an int* out-parameter.  They are
  ! passed as c_loc(...), so they too need the target attribute.
  integer(c_int), target :: major, minor, patch
  integer(c_int), target :: nc, ngroups, gsize

  integer(c_int), target              :: colors(n)
  integer(c_int), target, allocatable :: members(:)

  type(c_ptr)    :: result            ! opaque handle, owned by the library
  integer(c_int) :: ret
  integer        :: g, j, k

  ! -------------------------------------------------------------------------
  ! Version of SparseMatrixColorings.jl embedded in the library.  smc_version
  ! is the one void-returning entry point, hence a subroutine call; its three
  ! int* arguments are ordinary c_loc scalars.
  ! -------------------------------------------------------------------------
  call smc_version(c_loc(major), c_loc(minor), c_loc(patch))
  write(*,'(A,I0,A,I0,A,I0)') "SparseMatrixColorings ", major, ".", minor, ".", patch

  ! -------------------------------------------------------------------------
  ! Options.  smc_default_options() gives:
  !   structure     = SMC_NONSYMMETRIC
  !   partition     = SMC_COLUMN
  !   decompression = SMC_DIRECT
  !   order         = SMC_NATURAL
  !   index_base    = 0
  !   dtype         = SMC_FLOAT64
  ! Always start from that call, then override only the fields you care about —
  ! here the index base, so the CSC arrays above can be written the Fortran way.
  ! -------------------------------------------------------------------------
  opts = smc_default_options()
  opts%index_base = 1

  ! -------------------------------------------------------------------------
  ! Compute the coloring.  The result is an opaque handle owned by the library;
  ! it must be released with smc_result_free().
  !
  ! result_out is the one argument that is *not* passed by value: it is the
  ! void** out-parameter, declared type(c_ptr), intent(out), so `result` is
  ! given bare rather than through c_loc.  Everything else is c_loc(...).
  ! -------------------------------------------------------------------------
  result = c_null_ptr
  ret = smc_coloring(m, n, c_loc(colptr), c_loc(rowval), c_loc(opts), result)
  if (ret /= 0) then
    write(error_unit,'(A,I0,A)') "smc_coloring failed (", ret, ")"
    stop 1
  end if

  ! -------------------------------------------------------------------------
  ! Number of colors.
  ! -------------------------------------------------------------------------
  nc = 0
  ret = smc_ncolors(result, c_loc(nc))
  if (ret /= 0) then
    write(error_unit,'(A,I0,A)') "smc_ncolors failed (", ret, ")"
    ret = smc_result_free(result)
    stop 1
  end if
  write(*,'(A,I0)') "ncolors = ", nc

  ! -------------------------------------------------------------------------
  ! Per-column colors.  The buffer is caller-allocated and must hold at least n
  ! entries; the length is passed alongside it and checked before a single
  ! element is written, so a short buffer is refused with -3 rather than
  ! overrun.  Colors are integers in 1..ncolors; the value 0 is the neutral
  ! color, which only appears when opts%postprocessing is enabled and marks a
  ! column that needs no directional derivative at all.
  !
  ! Note that colors are labels, not indices: they are unaffected by
  ! opts%index_base, and so are indexed 1..n here purely because `colors` is a
  ! Fortran array.
  ! -------------------------------------------------------------------------
  ret = smc_column_colors(result, c_loc(colors), n)
  if (ret /= 0) then
    write(error_unit,'(A,I0,A)') "smc_column_colors failed (", ret, ")"
    ret = smc_result_free(result)
    stop 1
  end if

  write(*,'(A)', advance='no') "column colors = ["
  do j = 1, n
    write(*,'(A,I0)', advance='no') " ", colors(j)
  end do
  write(*,'(A)') " ]"

  ! -------------------------------------------------------------------------
  ! The same information seen the other way round: the groups of columns that
  ! share a color.  Group indices run over 1..ncolumn_groups and are always
  ! 1-based, independently of opts%index_base; the *members* are column indices
  ! written in the caller's index base, so with index_base = 1 they can be used
  ! to subscript a Fortran array directly.
  !
  ! Sizes are queried first so the caller can size its buffer exactly.
  ! -------------------------------------------------------------------------
  ngroups = 0
  ret = smc_ncolumn_groups(result, c_loc(ngroups))
  if (ret /= 0) then
    write(error_unit,'(A,I0,A)') "smc_ncolumn_groups failed (", ret, ")"
    ret = smc_result_free(result)
    stop 1
  end if

  do g = 1, ngroups
    gsize = 0
    ret = smc_column_group_size(result, int(g, c_int), c_loc(gsize))
    if (ret /= 0) then
      write(error_unit,'(A,I0,A,I0,A)') "smc_column_group_size failed for group ", g, &
                                        " (", ret, ")"
      ret = smc_result_free(result)
      stop 1
    end if

    ! c_loc requires a nonzero-sized object, so never allocate 0 elements even
    ! though a color class is always nonempty in practice.
    allocate(members(max(gsize, 1)))

    ret = smc_column_group(result, int(g, c_int), c_loc(members), gsize)
    if (ret /= 0) then
      write(error_unit,'(A,I0,A,I0,A)') "smc_column_group failed for group ", g, &
                                        " (", ret, ")"
      deallocate(members)
      ret = smc_result_free(result)
      stop 1
    end if

    write(*,'(A,I0,A)', advance='no') "group ", g, " = {"
    do k = 1, gsize
      write(*,'(A,I0)', advance='no') " ", members(k)
    end do
    write(*,'(A)') " }"

    deallocate(members)
  end do

  ! -------------------------------------------------------------------------
  ! Release the handle.  Using it after this point returns -4 (invalid or
  ! already-freed handle) rather than crashing, and freeing twice is safe.
  ! -------------------------------------------------------------------------
  ret = smc_result_free(result)
  if (ret /= 0) then
    write(error_unit,'(A,I0,A)') "smc_result_free failed (", ret, ")"
    stop 1
  end if

end program basic_coloring
