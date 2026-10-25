## Minimum vertex covers

#=
Vertex cover routines used by the post-processing to choose the colors kept for the trivial
two-colored structures (see `cover_trivial_structures!` in `postprocessing.jl`).

- `bipartite_vertex_cover!` computes a minimum vertex cover of a bipartite graph exactly, from a
  maximum matching and Kőnig's theorem.
- `greedy_cover!` computes a vertex cover of a general graph greedily, the minimum vertex cover
  problem being NP-hard in that case.
=#

# Compressed adjacency of the edges (a, b) from a to b, for a in 1:na
function _pair_adjacency(pairs::Vector{Tuple{Int,Int}}, na::Int)
    ptr = zeros(Int, na + 1)
    for (a, _) in pairs
        ptr[a + 1] += 1
    end
    ptr[1] = 1
    for a in 1:na
        ptr[a + 1] += ptr[a]
    end
    pos = ptr[1:na]
    adj = Vector{Int}(undef, length(pairs))
    for (a, b) in pairs
        adj[pos[a]] = b
        pos[a] += 1
    end
    return ptr, adj
end

# Same as `_pair_adjacency` for a bipartite graph with left vertices 1:nl and right vertices
# 1:nr, with duplicate edges removed by stamping the right vertices (no sort).
function _bipartite_adjacency(pairs::Vector{Tuple{Int,Int}}, nl::Int, nr::Int)
    ptr, adj = _pair_adjacency(pairs, nl)
    stamp = zeros(Int, nr)
    k = 0
    first = 1
    for a in 1:nl
        last = ptr[a + 1] - 1
        ptr[a] = k + 1
        for p in first:last
            b = adj[p]
            if stamp[b] != a
                stamp[b] = a
                k += 1
                adj[k] = b
            end
        end
        first = last + 1
    end
    ptr[nl + 1] = k + 1
    resize!(adj, k)
    return ptr, adj
end

## Maximum matching algorithms

"""
    AbstractMatchingAlgorithm

Algorithm computing a maximum matching of a bipartite graph, used by
`bipartite_vertex_cover!`. Both `KuhnMatching` and
`HopcroftKarpMatching` return a maximum matching; the minimum vertex cover derived from it
does not depend on the choice, only the running time does.
"""
abstract type AbstractMatchingAlgorithm end

"""
    KuhnMatching()

Kuhn's augmenting path algorithm: one search per left vertex, `O(|V| |E|)` time.
"""
struct KuhnMatching <: AbstractMatchingAlgorithm end

"""
    HopcroftKarpMatching()

Hopcroft-Karp algorithm: augments along a maximal set of shortest augmenting paths in each phase,
`O(|E| √|V|)` time.
"""
struct HopcroftKarpMatching <: AbstractMatchingAlgorithm end

## Minimum vertex cover of a bipartite graph

"""
    bipartite_vertex_cover!(cover_left, cover_right, ptr, adj, matching)

Set to `true` the entries of `cover_left` and `cover_right` that form a minimum vertex cover of
the bipartite graph with left vertices `1:length(ptr)-1`, right vertices `1:length(cover_right)`,
and neighbors `adj[ptr[a]:(ptr[a+1]-1)]` of the left vertex `a`.

The maximum matching is computed with `matching`, and the cover is read off it with Kőnig's
theorem: let `Z` be the set of vertices reachable from the unmatched left vertices by alternating
paths; then `(left \\ Z) ∪ (right ∩ Z)` is a minimum vertex cover. `Z` is the same for every
maximum matching, so the cover does not depend on `matching`.
"""
function bipartite_vertex_cover!(
    cover_left::AbstractVector{Bool},
    cover_right::AbstractVector{Bool},
    ptr::Vector{Int},
    adj::Vector{Int},
    matching::AbstractMatchingAlgorithm,
)
    nl = length(ptr) - 1
    nr = length(cover_right)
    match_left = zeros(Int, nl)
    match_right = zeros(Int, nr)
    maximum_matching!(match_left, match_right, ptr, adj, matching)

    # alternating breadth-first search from the unmatched left vertices
    in_z_left = zeros(Bool, nl)
    in_z_right = zeros(Bool, nr)
    queue = Vector{Int}(undef, nl)
    nq = 0
    for a in 1:nl
        if match_left[a] == 0 && ptr[a + 1] > ptr[a]
            in_z_left[a] = true
            nq += 1
            queue[nq] = a
        end
    end
    head = 1
    while head <= nq
        a = queue[head]
        head += 1
        for p in ptr[a]:(ptr[a + 1] - 1)
            b = adj[p]
            if !in_z_right[b] && match_left[a] != b
                in_z_right[b] = true
                a2 = match_right[b]  # b is matched, otherwise the matching is not maximum
                if !in_z_left[a2]
                    in_z_left[a2] = true
                    nq += 1
                    queue[nq] = a2
                end
            end
        end
    end

    for a in 1:nl
        if ptr[a + 1] > ptr[a] && !in_z_left[a]
            cover_left[a] = true
        end
    end
    for b in 1:nr
        if in_z_right[b]
            cover_right[b] = true
        end
    end
    return nothing
end

function maximum_matching!(match_left, match_right, ptr, adj, ::KuhnMatching)
    visited = zeros(Int, length(match_right))
    for a in 1:(length(ptr) - 1)
        _kuhn_augment!(a, a, visited, ptr, adj, match_left, match_right)
    end
    return nothing
end

# Look for an augmenting path from the left vertex a; `visited[b] == stamp` marks the right
# vertices already tried during the current search, which avoids resetting `visited`.
function _kuhn_augment!(a, stamp, visited, ptr, adj, match_left, match_right)
    for p in ptr[a]:(ptr[a + 1] - 1)
        b = adj[p]
        visited[b] == stamp && continue
        visited[b] = stamp
        if match_right[b] == 0 ||
            _kuhn_augment!(match_right[b], stamp, visited, ptr, adj, match_left, match_right)
            match_left[a] = b
            match_right[b] = a
            return true
        end
    end
    return false
end

function maximum_matching!(match_left, match_right, ptr, adj, ::HopcroftKarpMatching)
    nl = length(ptr) - 1
    dist = zeros(Int, nl)
    queue = Vector{Int}(undef, nl)
    while _hk_bfs!(dist, queue, ptr, adj, match_left, match_right)
        for a in 1:nl
            if match_left[a] == 0
                _hk_dfs!(a, dist, ptr, adj, match_left, match_right)
            end
        end
    end
    return nothing
end

# Layer the left vertices by alternating distance from the unmatched ones, and tell whether an
# augmenting path exists.
function _hk_bfs!(dist, queue, ptr, adj, match_left, match_right)
    nq = 0
    for a in eachindex(match_left)
        if match_left[a] == 0
            dist[a] = 0
            nq += 1
            queue[nq] = a
        else
            dist[a] = typemax(Int)
        end
    end
    found = false
    head = 1
    while head <= nq
        a = queue[head]
        head += 1
        for p in ptr[a]:(ptr[a + 1] - 1)
            a2 = match_right[adj[p]]
            if a2 == 0
                found = true
            elseif dist[a2] == typemax(Int)
                dist[a2] = dist[a] + 1
                nq += 1
                queue[nq] = a2
            end
        end
    end
    return found
end

# Augment along a shortest augmenting path from a, following the layers of `_hk_bfs!`.
function _hk_dfs!(a, dist, ptr, adj, match_left, match_right)
    for p in ptr[a]:(ptr[a + 1] - 1)
        b = adj[p]
        a2 = match_right[b]
        if a2 == 0 || (
            dist[a2] == dist[a] + 1 && _hk_dfs!(a2, dist, ptr, adj, match_left, match_right)
        )
            match_left[a] = b
            match_right[b] = a
            return true
        end
    end
    dist[a] = typemax(Int)  # dead end for this phase
    return false
end

## Greedy vertex cover of a general graph

# Greedy vertex cover of the undirected graph whose edges are `pairs`: repeatedly add the
# vertex incident to the largest number of uncovered edges, with degrees kept up to date.
function greedy_cover!(cover::Vector{Bool}, pairs::Vector{Tuple{Int,Int}})
    nv = length(cover)
    # undirected simple graph on the vertices, each pair stored in both directions
    sym_pairs = Tuple{Int,Int}[]
    for (a, b) in pairs
        push!(sym_pairs, (a, b))
        push!(sym_pairs, (b, a))
    end
    sym_pairs = unique!(sort!(sym_pairs))
    ptr, adj = _pair_adjacency(sym_pairs, nv)

    # dynamic degree = number of uncovered pairs incident to the vertex
    degree = [ptr[c + 1] - ptr[c] for c in 1:nv]
    dmax = maximum(degree; init=0)
    # buckets with lazy deletion: a vertex may appear in several buckets, only the one
    # matching its current degree is valid
    buckets = [Int[] for _ in 1:dmax]
    for c in nv:-1:1
        degree[c] > 0 && push!(buckets[degree[c]], c)
    end
    d = dmax
    while d > 0
        if isempty(buckets[d])
            d -= 1
            continue
        end
        c = pop!(buckets[d])
        (cover[c] || degree[c] != d) && continue
        # add the vertex c to the cover, which covers all its remaining pairs
        cover[c] = true
        degree[c] = 0
        for p in ptr[c]:(ptr[c + 1] - 1)
            c2 = adj[p]
            if !cover[c2]
                degree[c2] -= 1
                degree[c2] > 0 && push!(buckets[degree[c2]], c2)
            end
        end
    end
    return nothing
end
