# Guard for static compilation (`juliac --trim=safe`).
#
# The C and Fortran interfaces compile SMC.jl into a standalone shared library.
# The trimming verifier rejects any call whose result type is not fully inferred,
# so `coloring` must return a concrete type.

using SparseArrays
using SparseMatrixColorings
using Test

function _colored_result_type(problem, algo, ::Type{R}) where {R}
    return only(
        Base.return_types(
            (A, p, a, sp) -> coloring(A, p, a; decompression_eltype=R, symmetric_pattern=sp),
            (SparseMatrixCSC{Float64,Int}, typeof(problem), typeof(algo), Bool),
        ),
    )
end

const STATIC_COMBOS = [
    (ColoringProblem{:nonsymmetric,:column}(), GreedyColoringAlgorithm{:direct}()),
    (ColoringProblem{:nonsymmetric,:row}(), GreedyColoringAlgorithm{:direct}()),
    (ColoringProblem{:symmetric,:column}(), GreedyColoringAlgorithm{:direct}()),
    # The three below carry `decompression_eltype` in their result type, so they
    # are the ones that actually break without the annotation.
    (ColoringProblem{:symmetric,:column}(), GreedyColoringAlgorithm{:substitution}()),
    (ColoringProblem{:nonsymmetric,:bidirectional}(), GreedyColoringAlgorithm{:direct}()),
    (ColoringProblem{:nonsymmetric,:bidirectional}(), GreedyColoringAlgorithm{:substitution}()),
]

@testset "coloring infers a concrete result type" begin
    @testset "$(typeof(problem)) / $(typeof(algo)) / $R" for (problem, algo) in
                                                             STATIC_COMBOS,
        R in (Float32, Float64)

        RT = _colored_result_type(problem, algo, R)
        concrete = isconcretetype(RT)
        concrete || @error """
            `coloring` no longer infers a concrete result type.

            This breaks static compilation: `juliac --trim=safe` will reject the
            call with an "unresolved call" verifier error, and the C and Fortran
            interfaces will stop building.
            """ problem = typeof(problem) algorithm = typeof(algo) decompression_eltype = R inferred = RT
        @test concrete
    end
end
