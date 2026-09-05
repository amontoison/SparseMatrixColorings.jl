! compress_decompress.f90 — the full round trip: color, compress, decompress.
!
! This is the workflow a Jacobian/Hessian code actually uses:
!
!   1. color the columns of the sparsity pattern of A                 (smc_coloring)
!   2. evaluate one directional derivative per color; stacking them
!      side by side gives the compressed matrix B = A * S             (smc_compress)
!   3. scatter B back into the sparse structure of A                  (smc_decompress)
!
! Here step 2 is done by the library itself (smc_compress sums the columns of A
! that share a color), which lets us check that decompress(compress(A)) == A.
!
! A = tridiag(-1, 2, -1), 6x6, so three colors suffice and the compressed
! matrix is 6x3 instead of 6x6 — a 2x saving on derivative evaluations.
!
! This example is also the reference for the buffer-length discipline of the
! API.  Every numerical buffer is passed together with its length in elements,
! and every length can be queried from the handle before allocating:
!
!   nzval  : smc_nnz             -> nnz            elements
!   Br, Bc : smc_compressed_size -> Br_rows*Br_cols, Bc_rows*Bc_cols elements
!   A_out  : smc_size            -> m*n            elements
!
! A buffer shorter than its required minimum is rejected with -3 before a
! single element is read or written; the last section below demonstrates it.
!
! Dense matrices crossing the interface are column-major, which is Fortran's
! own layout — so Bc and A below are plain 2D arrays, handed over with a single
! c_loc and indexed Bc(i,j) / A(i,j) with no manual i + j*rows arithmetic.
!
! Like basic_coloring.f90 this sets opts%index_base = 1, so colptr and rowval
! are 1-based and rowval(k) subscripts A directly.
!
! Compile (after building libsmc with juliac — see interfaces/README.md):
!
!   gfortran -O2 -o interfaces/build/compress_decompress_f \
!       interfaces/examples/Fortran/compress_decompress.f90 \
!       -I interfaces/include \
!       interfaces/build/lib/libsmc.so \
!       -Wl,-rpath,'$ORIGIN/lib'
!
! -I points gfortran at the directory holding smc.f90.  The rpath assumes the
! binary sits in interfaces/build/, next to lib/; on macOS use
! -Wl,-rpath,@loader_path/lib instead.
!
! Expected output:
!   ncolors = 3
!   pattern: 6 x 6, 16 stored entries
!   compressed matrix B (6 x 3):
!      2.00  -1.00   0.00
!     -1.00   2.00  -1.00
!     -1.00  -1.00   2.00
!      2.00  -1.00  -1.00
!     -1.00   2.00  -1.00
!      0.00  -1.00   2.00
!   max error on the 16 stored entries: 0.000E+00
!   max leakage outside the pattern:    0.000E+00
!   short nzval -> -3, short Bc -> -3, short A_out -> -3
!
! Exit code: 0 on success, 1 if any call fails or the round trip is inexact.

program compress_decompress
  use iso_c_binding
  use iso_fortran_env, only: error_unit
  implicit none
  include 'smc.f90'          ! <- after implicit none

  ! -------------------------------------------------------------------------
  ! Problem data: A = tridiag(-1, 2, -1), 6x6, compressed sparse column,
  ! 1-based to match opts%index_base = 1 below.
  !
  ! nzval follows exactly the same ordering as rowval: nzval(k) is the value of
  ! the entry stored at row rowval(k).  Coloring ignores it entirely; only
  ! smc_compress reads it.
  !
  ! All four arrays carry target, since their addresses are taken with c_loc.
  ! -------------------------------------------------------------------------

  integer(c_int), parameter :: M_ROWS = 6   ! number of rows
  integer(c_int), parameter :: N_COLS = 6   ! number of columns
  integer(c_int), parameter :: NNZ    = 16  ! number of stored entries

  integer(c_int), target :: colptr(N_COLS+1) = [ 1, 3, 6, 9, 12, 15, 17 ]
  integer(c_int), target :: rowval(NNZ)      = [ 1, 2,        &
                                                 1, 2, 3,     &
                                                 2, 3, 4,     &
                                                 3, 4, 5,     &
                                                 4, 5, 6,     &
                                                 5, 6 ]
  real(c_double), target :: nzval(NNZ)       = [  2.0_c_double, -1.0_c_double,                 &
                                                 -1.0_c_double,  2.0_c_double, -1.0_c_double,  &
                                                 -1.0_c_double,  2.0_c_double, -1.0_c_double,  &
                                                 -1.0_c_double,  2.0_c_double, -1.0_c_double,  &
                                                 -1.0_c_double,  2.0_c_double, -1.0_c_double,  &
                                                 -1.0_c_double,  2.0_c_double ]

  type(SmcColoringOptions), target :: opts

  ! Scalars filled by the library through int* out-parameters.
  integer(c_int), target :: nc, nnz_q, m, n
  integer(c_int), target :: Br_rows, Br_cols, Bc_rows, Bc_cols

  ! Dense buffers, allocated from the queried sizes and never from a guess.
  real(c_double), target, allocatable :: Bc(:,:)   ! compressed   (Bc_rows x Bc_cols)
  real(c_double), target, allocatable :: A(:,:)    ! decompressed (m x n)
  logical,                allocatable :: in_pattern(:,:)

  ! Lengths are element counts, never byte counts, and use c_size_t because
  ! m*n overflows a 32-bit int for perfectly ordinary dimensions.
  integer(c_size_t) :: nzval_len, Bc_len, A_len

  type(c_ptr)    :: result
  integer(c_int) :: ret, short_nzval, short_Bc, short_A
  integer        :: i, j, k
  real(c_double) :: max_err, max_leak, err, leak

  ! -------------------------------------------------------------------------
  ! Defaults: nonsymmetric / column / direct / natural, Float64; we switch the
  ! index base to 1 for Fortran.
  !
  ! dtype is what fixes the element type of the compressed and decompressed
  ! buffers below: SMC_FLOAT64 means real(c_double), SMC_FLOAT32 would mean
  ! real(c_float).  It also fixes the unit of every *_len argument, which
  ! counts elements of that type — never bytes.
  ! -------------------------------------------------------------------------
  opts = smc_default_options()
  opts%index_base = 1

  result = c_null_ptr
  ret = smc_coloring(M_ROWS, N_COLS, c_loc(colptr), c_loc(rowval), c_loc(opts), result)
  if (ret /= 0) then
    write(error_unit,'(A,I0,A)') "smc_coloring failed (", ret, ")"
    stop 1
  end if

  nc = 0
  ret = smc_ncolors(result, c_loc(nc))
  if (ret /= 0) then
    write(error_unit,'(A,I0,A)') "smc_ncolors failed (", ret, ")"
    ret = smc_result_free(result)
    stop 1
  end if
  write(*,'(A,I0)') "ncolors = ", nc

  ! -------------------------------------------------------------------------
  ! Ask the handle for every size we are about to allocate.  A caller that
  ! received the handle from elsewhere and has none of the CSC arrays at hand
  ! can still do this: the result remembers the pattern it was built from.
  !
  !   smc_nnz  : length nzval must have in smc_compress
  !   smc_size : the m and n of that pattern, so A must hold m*n elements
  ! -------------------------------------------------------------------------
  nnz_q = 0
  ret = smc_nnz(result, c_loc(nnz_q))
  if (ret /= 0) then
    write(error_unit,'(A,I0,A)') "smc_nnz failed (", ret, ")"
    ret = smc_result_free(result)
    stop 1
  end if

  m = 0
  n = 0
  ret = smc_size(result, c_loc(m), c_loc(n))
  if (ret /= 0) then
    write(error_unit,'(A,I0,A)') "smc_size failed (", ret, ")"
    ret = smc_result_free(result)
    stop 1
  end if
  write(*,'(A,I0,A,I0,A,I0,A)') "pattern: ", m, " x ", n, ", ", nnz_q, " stored entries"

  ! Our static data must match what the handle reports, or the lengths we are
  ! about to pass would be a lie.
  if (m /= M_ROWS .or. n /= N_COLS .or. nnz_q /= NNZ) then
    write(error_unit,'(A)') "handle describes a different pattern"
    ret = smc_result_free(result)
    stop 1
  end if

  ! -------------------------------------------------------------------------
  ! Shape of the compressed matrix.
  !
  ! With a bidirectional partition the compression is split in two blocks, a
  ! row block Br and a column block Bc.  For the column partition used here
  ! only Bc is meaningful and the Br dimensions come back as 0, so we pass
  ! c_null_ptr wherever Br is expected — with a length of 0 to match.
  ! -------------------------------------------------------------------------
  Br_rows = 0
  Br_cols = 0
  Bc_rows = 0
  Bc_cols = 0
  ret = smc_compressed_size(result, c_loc(Br_rows), c_loc(Br_cols), &
                                    c_loc(Bc_rows), c_loc(Bc_cols))
  if (ret /= 0) then
    write(error_unit,'(A,I0,A)') "smc_compressed_size failed (", ret, ")"
    ret = smc_result_free(result)
    stop 1
  end if

  ! Allocate exactly what the queries reported — nothing is hard-coded here.
  ! Bc is m-by-ncolors for a column partition.  Because Fortran stores arrays
  ! column-major, a 2D array is already in the layout the library expects.
  allocate(Bc(Bc_rows, Bc_cols))
  allocate(A(m, n))
  Bc = 0.0_c_double
  A  = 0.0_c_double

  nzval_len = int(nnz_q,   c_size_t)
  Bc_len    = int(Bc_rows, c_size_t) * int(Bc_cols, c_size_t)
  A_len     = int(m,       c_size_t) * int(n,       c_size_t)

  ! -------------------------------------------------------------------------
  ! Compress: B(:, c) = sum of the columns of A colored c.
  ! In a real AD code this buffer would instead be filled by evaluating one
  ! directional derivative per color.
  !
  ! Each buffer is followed by its length in elements.  Br is unused by a
  ! column partition, so it is c_null_ptr with length 0.
  ! -------------------------------------------------------------------------
  ret = smc_compress(result,                        &
                     c_loc(nzval), nzval_len,       &
                     c_null_ptr,   0_c_size_t,      &  ! Br unused by a column partition
                     c_loc(Bc),    Bc_len)
  if (ret /= 0) then
    write(error_unit,'(A,I0,A)') "smc_compress failed (", ret, ")"
    ret = smc_result_free(result)
    stop 1
  end if

  write(*,'(A,I0,A,I0,A)') "compressed matrix B (", Bc_rows, " x ", Bc_cols, "):"
  do i = 1, Bc_rows
    do j = 1, Bc_cols
      write(*,'(F7.2)', advance='no') Bc(i,j)
    end do
    write(*,*)
  end do

  ! -------------------------------------------------------------------------
  ! Decompress: recover the full m-by-n dense matrix A from B.  Again the
  ! column-major buffer is just a Fortran 2D array, so A(i,j) is the entry the
  ! C example has to spell Ad[i + j*m].
  ! -------------------------------------------------------------------------
  ret = smc_decompress(result,                    &
                       c_null_ptr, 0_c_size_t,    &  ! Br unused by a column partition
                       c_loc(Bc),  Bc_len,        &
                       c_loc(A),   A_len)
  if (ret /= 0) then
    write(error_unit,'(A,I0,A)') "smc_decompress failed (", ret, ")"
    ret = smc_result_free(result)
    stop 1
  end if

  ! -------------------------------------------------------------------------
  ! Check 1: every stored entry of A came back exactly.
  ! Check 2: nothing leaked into positions outside the sparsity pattern.
  !
  ! With index_base = 1 the CSC arrays are already Fortran indices, so the
  ! traversal below needs no shifting at all: column j runs over the entries
  ! colptr(j) .. colptr(j+1)-1, and rowval(k) subscripts A directly.
  ! -------------------------------------------------------------------------
  allocate(in_pattern(m, n))
  in_pattern = .false.
  max_err = 0.0_c_double

  do j = 1, N_COLS
    do k = colptr(j), colptr(j+1) - 1
      i = rowval(k)
      in_pattern(i, j) = .true.
      err = abs(A(i,j) - nzval(k))
      if (err > max_err) max_err = err
    end do
  end do
  write(*,'(A,I0,A,ES9.3)') "max error on the ", NNZ, " stored entries: ", max_err

  max_leak = 0.0_c_double
  do j = 1, N_COLS
    do i = 1, M_ROWS
      if (.not. in_pattern(i,j)) then
        leak = abs(A(i,j))
        if (leak > max_leak) max_leak = leak
      end if
    end do
  end do
  write(*,'(A,ES9.3)') "max leakage outside the pattern:    ", max_leak

  ! -------------------------------------------------------------------------
  ! The lengths are not decoration.  Understating any of them is refused with
  ! -3, and the library reads and writes nothing before it refuses — which is
  ! exactly the guarantee that lets a caller reuse a handle safely.
  ! -------------------------------------------------------------------------
  short_nzval = smc_compress(result, c_loc(nzval), nzval_len - 1, &
                             c_null_ptr, 0_c_size_t, c_loc(Bc), Bc_len)
  short_Bc    = smc_compress(result, c_loc(nzval), nzval_len,     &
                             c_null_ptr, 0_c_size_t, c_loc(Bc), Bc_len - 1)
  short_A     = smc_decompress(result, c_null_ptr, 0_c_size_t,    &
                               c_loc(Bc), Bc_len, c_loc(A), A_len - 1)
  write(*,'(A,I0,A,I0,A,I0)') "short nzval -> ", short_nzval, &
                              ", short Bc -> ",  short_Bc,    &
                              ", short A_out -> ", short_A

  ! -------------------------------------------------------------------------
  ! Release the handle; the Fortran arrays are deallocated on the way out.
  ! -------------------------------------------------------------------------
  ret = smc_result_free(result)
  if (ret /= 0) then
    write(error_unit,'(A,I0,A)') "smc_result_free failed (", ret, ")"
    stop 1
  end if

  deallocate(in_pattern, Bc, A)

  ! Direct decompression is an exact scatter, so both errors must be zero.
  ! (Both are absolute values, hence nonnegative: > 0 is the exact-equality
  ! test, spelled so that no compiler warns about comparing reals.)
  if (max_err > 0.0_c_double .or. max_leak > 0.0_c_double) then
    write(error_unit,'(A)') "round trip did not reproduce A"
    stop 1
  end if

  if (short_nzval /= -3 .or. short_Bc /= -3 .or. short_A /= -3) then
    write(error_unit,'(A)') "a short buffer was not rejected"
    stop 1
  end if

end program compress_decompress
