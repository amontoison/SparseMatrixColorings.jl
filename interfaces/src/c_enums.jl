# c_enums.jl — Julia mirror of the C types declared in include/smc.h.
#
# The enumerator *values* live in scripts/coloring_table.jl (the single source
# of truth shared with the header generator) and are already in scope here:
#
#   SmcDataType       SMC_FLOAT64 = 0, SMC_FLOAT32 = 1
#   SmcStructure      SMC_NONSYMMETRIC = 0, SMC_SYMMETRIC = 1
#   SmcPartition      SMC_COLUMN = 0, SMC_ROW = 1, SMC_BIDIRECTIONAL = 2
#   SmcDecompression  SMC_DIRECT = 0, SMC_SUBSTITUTION = 1
#   SmcOrder          SMC_NATURAL = 0 ... SMC_DYNAMIC_LARGEST_FIRST = 4
#
# The GC roots for the opaque result handles are the typed stores of
# c_stores.jl.

# ---------------------------------------------------------------------------
# SmcColoringOptionsC — Julia mirror of the SmcColoringOptions C struct.
#
# Eight `Cint` fields in exactly the order of DESIGN.md section 2 and of
# include/smc.h.  All fields are 4 bytes, so the C and Julia layouts agree with
# no padding on every supported platform, and the struct is `isbits`: it can be
# `unsafe_load`ed straight from the caller's pointer and returned by value from
# `smc_default_options`.
#
# Changing the order of these fields is an ABI break.
# ---------------------------------------------------------------------------
struct SmcColoringOptionsC
  structure         :: Cint   # SmcStructure      — default SMC_NONSYMMETRIC
  partition         :: Cint   # SmcPartition      — default SMC_COLUMN
  decompression     :: Cint   # SmcDecompression  — default SMC_DIRECT
  order             :: Cint   # SmcOrder          — default SMC_NATURAL
  postprocessing    :: Cint   # 0/1 — neutral color 0 where possible
  symmetric_pattern :: Cint   # 0/1 — assert the sparsity pattern is symmetric
  index_base        :: Cint   # 0 or 1 — index base of colptr / rowval / groups
  dtype             :: Cint   # SmcDataType — element type of compress/decompress
end

# The defaults documented in DESIGN.md section 2 and returned both by
# `smc_default_options` and whenever a NULL options pointer is passed.
const SMC_DEFAULT_OPTIONS = SmcColoringOptionsC(
  SMC_NONSYMMETRIC,   # structure
  SMC_COLUMN,         # partition
  SMC_DIRECT,         # decompression
  SMC_NATURAL,        # order
  Cint(0),            # postprocessing
  Cint(0),            # symmetric_pattern
  Cint(0),            # index_base (0 = C-natural)
  SMC_FLOAT64,        # dtype
)

# Every field of an options struct is in range.  Anything out of range is an
# invalid argument (-3), which is a different failure from a well-formed but
# unsupported combination (-2).  `postprocessing` and `symmetric_pattern` are
# read as booleans, so any nonzero value is accepted for them.
function _valid_options(o::SmcColoringOptionsC)
  Cint(0) <= o.structure <= SMC_MAX_STRUCTURE || return false
  Cint(0) <= o.partition <= SMC_MAX_PARTITION || return false
  Cint(0) <= o.decompression <= SMC_MAX_DECOMPRESSION || return false
  Cint(0) <= o.order <= SMC_MAX_ORDER || return false
  Cint(0) <= o.dtype <= SMC_MAX_DTYPE || return false
  o.index_base == Cint(0) || o.index_base == Cint(1) || return false
  return true
end
