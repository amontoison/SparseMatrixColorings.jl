# coloring_table.jl — single source of truth for the libsmc combination tables.
#
# This file is pure data: no `using`, no dependency on SparseMatrixColorings, no
# reference to any type of the package.  That is what lets it be loaded by two
# very different consumers:
#
#   * `interfaces/src/c_stores.jl` (through `interfaces/src/LibSMC.jl`), which
#     needs the enum values and the combo keys that index the typed handle
#     stores;
#   * `interfaces/scripts/generate_header.jl`, which loads it into a bare
#     `Module` — with no package environment at all — to cross-check the
#     enumerators it writes into `include/smc.h`.
#
# Everything here is a verbatim transcription of DESIGN.md sections 2 and 3.
# The *values* are part of the ABI: entry `i` of each enum table has the value
# `i - 1`, and those numbers must never move.

# ---------------------------------------------------------------------------
# Enumerators (DESIGN.md section 2).  Each entry is
#   (C enumerator, Julia counterpart as it appears in LibSMC.jl)
# and its value is its 0-based position in the table.
# ---------------------------------------------------------------------------

const DTYPES = [
  ("SMC_FLOAT64", "Float64"),
  ("SMC_FLOAT32", "Float32"),
]

const STRUCTURES = [
  ("SMC_NONSYMMETRIC", ":nonsymmetric"),
  ("SMC_SYMMETRIC", ":symmetric"),
]

const PARTITIONS = [
  ("SMC_COLUMN", ":column"),
  ("SMC_ROW", ":row"),
  ("SMC_BIDIRECTIONAL", ":bidirectional"),
]

const DECOMPRESSIONS = [
  ("SMC_DIRECT", ":direct"),
  ("SMC_SUBSTITUTION", ":substitution"),
]

# `RandomOrder` is deliberately excluded from v1: it carries an `AbstractRNG`,
# which is untested under `--trim=safe`.  Every order below is a singleton, so
# each entry can be written as a literal constructor inside an `if` branch
# (constraint C1 of DESIGN.md section 0).
const ORDERS = [
  ("SMC_NATURAL", "NaturalOrder()"),
  ("SMC_LARGEST_FIRST", "LargestFirst()"),
  ("SMC_SMALLEST_LAST", "DynamicDegreeBasedOrder{:back,:high2low,false}()"),
  ("SMC_INCIDENCE_DEGREE", "DynamicDegreeBasedOrder{:back,:low2high,false}()"),
  ("SMC_DYNAMIC_LARGEST_FIRST", "DynamicDegreeBasedOrder{:forward,:low2high,false}()"),
]

# ---------------------------------------------------------------------------
# Enum values, as `Cint` so they can be compared with the fields of
# `SmcColoringOptions` without any promotion.
# ---------------------------------------------------------------------------

const SMC_FLOAT64 = Cint(0)
const SMC_FLOAT32 = Cint(1)

const SMC_NONSYMMETRIC = Cint(0)
const SMC_SYMMETRIC = Cint(1)

const SMC_COLUMN = Cint(0)
const SMC_ROW = Cint(1)
const SMC_BIDIRECTIONAL = Cint(2)

const SMC_DIRECT = Cint(0)
const SMC_SUBSTITUTION = Cint(1)

const SMC_NATURAL = Cint(0)
const SMC_LARGEST_FIRST = Cint(1)
const SMC_SMALLEST_LAST = Cint(2)
const SMC_INCIDENCE_DEGREE = Cint(3)
const SMC_DYNAMIC_LARGEST_FIRST = Cint(4)

# Largest admissible value of each option, used by the -3 range checks.
const SMC_MAX_DTYPE = SMC_FLOAT32
const SMC_MAX_STRUCTURE = SMC_SYMMETRIC
const SMC_MAX_PARTITION = SMC_BIDIRECTIONAL
const SMC_MAX_DECOMPRESSION = SMC_SUBSTITUTION
const SMC_MAX_ORDER = SMC_DYNAMIC_LARGEST_FIRST

# ---------------------------------------------------------------------------
# Combo key (DESIGN.md section 3)
#
#   key = structure * 16 + partition * 4 + decompression * 2 + dtype
#
# The key fits in a `UInt8` (largest supported value is 0x13) and is what the
# handle -> key store records, so that every later call can recover the
# concrete type of the result it was given.
# ---------------------------------------------------------------------------

combo_key(structure, partition, decompression, dtype) =
  UInt8(structure * 16 + partition * 4 + decompression * 2 + dtype)

# Key that no combination can produce: "this handle is unknown" (-> return -4).
const KEY_INVALID = 0xff

const KEY_NS_COL_DIRECT_F64 = combo_key(SMC_NONSYMMETRIC, SMC_COLUMN, SMC_DIRECT, SMC_FLOAT64)
const KEY_NS_COL_DIRECT_F32 = combo_key(SMC_NONSYMMETRIC, SMC_COLUMN, SMC_DIRECT, SMC_FLOAT32)
const KEY_NS_ROW_DIRECT_F64 = combo_key(SMC_NONSYMMETRIC, SMC_ROW, SMC_DIRECT, SMC_FLOAT64)
const KEY_NS_ROW_DIRECT_F32 = combo_key(SMC_NONSYMMETRIC, SMC_ROW, SMC_DIRECT, SMC_FLOAT32)
const KEY_NS_BID_DIRECT_F64 = combo_key(SMC_NONSYMMETRIC, SMC_BIDIRECTIONAL, SMC_DIRECT, SMC_FLOAT64)
const KEY_NS_BID_DIRECT_F32 = combo_key(SMC_NONSYMMETRIC, SMC_BIDIRECTIONAL, SMC_DIRECT, SMC_FLOAT32)
const KEY_NS_BID_SUBST_F64 = combo_key(SMC_NONSYMMETRIC, SMC_BIDIRECTIONAL, SMC_SUBSTITUTION, SMC_FLOAT64)
const KEY_NS_BID_SUBST_F32 = combo_key(SMC_NONSYMMETRIC, SMC_BIDIRECTIONAL, SMC_SUBSTITUTION, SMC_FLOAT32)
const KEY_SYM_COL_DIRECT_F64 = combo_key(SMC_SYMMETRIC, SMC_COLUMN, SMC_DIRECT, SMC_FLOAT64)
const KEY_SYM_COL_DIRECT_F32 = combo_key(SMC_SYMMETRIC, SMC_COLUMN, SMC_DIRECT, SMC_FLOAT32)
const KEY_SYM_COL_SUBST_F64 = combo_key(SMC_SYMMETRIC, SMC_COLUMN, SMC_SUBSTITUTION, SMC_FLOAT64)
const KEY_SYM_COL_SUBST_F32 = combo_key(SMC_SYMMETRIC, SMC_COLUMN, SMC_SUBSTITUTION, SMC_FLOAT32)

# The same key with the `dtype` bit cleared: the three (structure, partition,
# decompression) triples that `smc_fast_coloring` needs, since it never builds
# a result object and therefore never looks at `dtype`.
const BKEY_NS_COL_DIRECT = KEY_NS_COL_DIRECT_F64
const BKEY_NS_ROW_DIRECT = KEY_NS_ROW_DIRECT_F64
const BKEY_NS_BID_DIRECT = KEY_NS_BID_DIRECT_F64
const BKEY_NS_BID_SUBST = KEY_NS_BID_SUBST_F64
const BKEY_SYM_COL_DIRECT = KEY_SYM_COL_DIRECT_F64
const BKEY_SYM_COL_SUBST = KEY_SYM_COL_SUBST_F64

# ---------------------------------------------------------------------------
# The nine typed stores of DESIGN.md section 3.
#
# Each entry is
#   (store number, store name, [keys routed to it], structure, partition,
#    decompression, dtype or "both", concrete value type as measured by
#    running Julia on a `SparseMatrixCSC{Float64,Int64}` input)
#
# `decompression_eltype` does not appear in the Column / Row / StarSet result
# types, so those three are dtype-independent: nine stores serve twelve keys.
# `interfaces/src/c_stores.jl` spells the types out again in Julia syntax; this
# table is the reference they are checked against.
#
# `SubArray{...}` below is exactly
#   SubArray{Int64,1,Vector{Int64},Tuple{UnitRange{Int64}},true}
# ---------------------------------------------------------------------------

const COMBOS = [
  (1, "store_ns_col_direct",
   [KEY_NS_COL_DIRECT_F64, KEY_NS_COL_DIRECT_F32],
   "SMC_NONSYMMETRIC", "SMC_COLUMN", "SMC_DIRECT", "both",
   "ColumnColoringResult{SparseMatrixCSC{Float64,Int64}, Int64, BipartiteGraph{Int64}, Vector{Int64}, Vector{SubArray{...}}, Vector{Int64}, Nothing}"),
  (2, "store_ns_row_direct",
   [KEY_NS_ROW_DIRECT_F64, KEY_NS_ROW_DIRECT_F32],
   "SMC_NONSYMMETRIC", "SMC_ROW", "SMC_DIRECT", "both",
   "RowColoringResult{SparseMatrixCSC{Float64,Int64}, Int64, BipartiteGraph{Int64}, Vector{Int64}, Vector{SubArray{...}}, Vector{Int64}, Nothing}"),
  (3, "store_sym_col_direct",
   [KEY_SYM_COL_DIRECT_F64, KEY_SYM_COL_DIRECT_F32],
   "SMC_SYMMETRIC", "SMC_COLUMN", "SMC_DIRECT", "both",
   "StarSetColoringResult{SparseMatrixCSC{Float64,Int64}, Int64, AdjacencyGraph{Int64,false}, Vector{Int64}, Vector{SubArray{...}}, Vector{Int64}, Nothing}"),
  (4, "store_sym_col_subst_f64",
   [KEY_SYM_COL_SUBST_F64],
   "SMC_SYMMETRIC", "SMC_COLUMN", "SMC_SUBSTITUTION", "SMC_FLOAT64",
   "TreeSetColoringResult{SparseMatrixCSC{Float64,Int64}, Int64, AdjacencyGraph{Int64,false}, Vector{SubArray{...}}, Float64}"),
  (5, "store_sym_col_subst_f32",
   [KEY_SYM_COL_SUBST_F32],
   "SMC_SYMMETRIC", "SMC_COLUMN", "SMC_SUBSTITUTION", "SMC_FLOAT32",
   "TreeSetColoringResult{SparseMatrixCSC{Float64,Int64}, Int64, AdjacencyGraph{Int64,false}, Vector{SubArray{...}}, Float32}"),
  (6, "store_ns_bid_direct_f64",
   [KEY_NS_BID_DIRECT_F64],
   "SMC_NONSYMMETRIC", "SMC_BIDIRECTIONAL", "SMC_DIRECT", "SMC_FLOAT64",
   "BicoloringResult{SparseMatrixCSC{Float64,Int64}, Int64, AdjacencyGraph{Int64,true}, :direct, Vector{SubArray{...}}, StarSetColoringResult{SparsityPatternCSC{Int64}, Int64, AdjacencyGraph{Int64,true}, Vector{Int64}, Vector{SubArray{...}}, Vector{Int64}, Nothing}, Float64}"),
  (7, "store_ns_bid_direct_f32",
   [KEY_NS_BID_DIRECT_F32],
   "SMC_NONSYMMETRIC", "SMC_BIDIRECTIONAL", "SMC_DIRECT", "SMC_FLOAT32",
   "BicoloringResult{SparseMatrixCSC{Float64,Int64}, Int64, AdjacencyGraph{Int64,true}, :direct, Vector{SubArray{...}}, StarSetColoringResult{SparsityPatternCSC{Int64}, Int64, AdjacencyGraph{Int64,true}, Vector{Int64}, Vector{SubArray{...}}, Vector{Int64}, Nothing}, Float32}"),
  (8, "store_ns_bid_subst_f64",
   [KEY_NS_BID_SUBST_F64],
   "SMC_NONSYMMETRIC", "SMC_BIDIRECTIONAL", "SMC_SUBSTITUTION", "SMC_FLOAT64",
   "BicoloringResult{SparseMatrixCSC{Float64,Int64}, Int64, AdjacencyGraph{Int64,true}, :substitution, Vector{SubArray{...}}, TreeSetColoringResult{SparsityPatternCSC{Int64}, Int64, AdjacencyGraph{Int64,true}, Vector{SubArray{...}}, Float64}, Float64}"),
  (9, "store_ns_bid_subst_f32",
   [KEY_NS_BID_SUBST_F32],
   "SMC_NONSYMMETRIC", "SMC_BIDIRECTIONAL", "SMC_SUBSTITUTION", "SMC_FLOAT32",
   "BicoloringResult{SparseMatrixCSC{Float64,Int64}, Int64, AdjacencyGraph{Int64,true}, :substitution, Vector{SubArray{...}}, TreeSetColoringResult{SparsityPatternCSC{Int64}, Int64, AdjacencyGraph{Int64,true}, Vector{SubArray{...}}, Float32}, Float32}"),
]

# Every other combination is rejected with -2.  In particular
# `(nonsymmetric, column|row, substitution)` is not a `ColoringProblem` that
# SparseMatrixColorings supports, and neither is any symmetric row or
# bidirectional partition.
const SUPPORTED_KEYS = [k for combo in COMBOS for k in combo[3]]
