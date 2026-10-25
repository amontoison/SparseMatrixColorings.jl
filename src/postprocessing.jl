## Postprocessing

function postprocess!(
    color::AbstractVector{<:Integer},
    star_or_tree_set::Union{StarSet,TreeSet},
    g::AdjacencyGraph,
    offsets::AbstractVector{<:Integer},
    postprocessing_minimizes::Symbol,
)
    # flag which colors are actually used during decompression
    nb_colors = maximum(color; init=0)
    color_used = zeros(Bool, nb_colors)

    # nonzero diagonal coefficients force the use of their respective color (there can be no neutral colors if the diagonal is fully nonzero)
    bicoloring = augmented_graph(g)
    if !bicoloring
        S = pattern(g)
        for i in axes(S, 1)
            if !iszero(S[i, i])
                color_used[color[i]] = true
            end
        end
    end

    # When bicoloring is false, column_color_used and row_color_used point to the same memory
    row_color_used = bicoloring ? zeros(Bool, nb_colors) : color_used
    column_color_used = color_used

    if star_or_tree_set isa StarSet
        # star_or_tree_set is a StarSet
        postprocess_with_star_set!(
            bicoloring,
            g,
            row_color_used,
            column_color_used,
            color,
            star_or_tree_set,
            postprocessing_minimizes,
        )
    else
        # star_or_tree_set is a TreeSet
        postprocess_with_tree_set!(
            bicoloring,
            row_color_used,
            column_color_used,
            color,
            star_or_tree_set,
            postprocessing_minimizes,
        )
    end

    # if at least one of the colors is not used, modify the color assignments of vertices
    has_neutral_color = if bicoloring
        any(!, row_color_used) || any(!, column_color_used)
    else
        any(!, color_used)
    end

    if has_neutral_color
        # size of the original matrix on which we want to perform coloring or bicoloring
        (m, n) = g.original_size

        # count the number of unused colors
        num_unused_colors = 0

        # count how many color indices are skipped before each color,
        # in order to compact the color indexing after removing unused colors
        for ci in 1:nb_colors
            ci_required =
                bicoloring ? row_color_used[ci] || column_color_used[ci] : color_used[ci]
            if ci_required
                offsets[ci] = num_unused_colors
            else
                num_unused_colors += 1
            end
        end

        # replace unused colors by the neutral color and compact the remaining color indices
        for i in eachindex(color)
            ci = color[i]
            ci_used = (i ≤ n) ? column_color_used[ci] : row_color_used[ci]
            if !ci_used
                # assign the neutral color
                color[i] = 0
            else
                # remap the color to not have any gap
                color[i] -= offsets[ci]
            end
        end
    end
    return color
end

function postprocess_with_star_set!(
    bicoloring::Bool,
    g::AdjacencyGraph,
    row_color_used::Vector{Bool},
    column_color_used::Vector{Bool},
    color::AbstractVector{<:Integer},
    star_set::StarSet,
    postprocessing_minimizes::Symbol=:all_colors,
)
    S = pattern(g)
    edge_to_index = edge_indices(g)

    # only the colors of the hubs are used
    (; star, hub) = star_set
    nb_trivial_stars = 0

    # size of the original matrix on which we want to perform bicoloring
    (m, n) = g.original_size

    # Iterate through all non-trivial stars
    for s in eachindex(hub)
        h = hub[s]
        if h > 0
            if h ≤ n
                column_color_used[color[h]] = true
            else
                row_color_used[color[h]] = true
            end
        else
            nb_trivial_stars += 1
        end
    end

    # Process the trivial stars (if any)
    if nb_trivial_stars > 0
        # Color pairs (color of i, color of j) of the trivial stars {i, j} with i > j whose
        # hub is not determined yet: neither color is already required on its side.
        undecided_pairs = Tuple{Int,Int}[]

        rvS = rowvals(S)
        for j in axes(S, 2)
            for k in nzrange(S, j)
                i = rvS[k]
                if i > j
                    s = star[edge_to_index[k]]
                    if hub[s] < 0
                        if row_color_used[color[i]]
                            # The vertex i is already a hub in a non-trivial star
                            hub[s] = i
                        elseif column_color_used[color[j]]
                            # The vertex j is already a hub in a non-trivial star
                            hub[s] = j
                        else
                            push!(undecided_pairs, (color[i], color[j]))
                        end
                    end
                end
            end
        end

        if !isempty(undecided_pairs)
            # Decide which colors to keep for the remaining trivial stars
            cover_trivial_structures!(
                row_color_used,
                column_color_used,
                undecided_pairs,
                bicoloring,
                postprocessing_minimizes,
            )

            # Every remaining trivial star now has an endpoint with a used color: make it the hub
            for j in axes(S, 2)
                for k in nzrange(S, j)
                    i = rvS[k]
                    if i > j
                        s = star[edge_to_index[k]]
                        if hub[s] < 0
                            hub[s] = row_color_used[color[i]] ? i : j
                        end
                    end
                end
            end
        end
    end
    return nothing
end

function postprocess_with_tree_set!(
    bicoloring::Bool,
    row_color_used::Vector{Bool},
    column_color_used::Vector{Bool},
    color::AbstractVector{<:Integer},
    tree_set::TreeSet,
    postprocessing_minimizes::Symbol=:all_colors,
)
    # only the colors of non-leaf vertices are used
    (; reverse_bfs_orders, is_star, tree_edge_indices, nt) = tree_set
    nb_trivial_trees = 0

    # Iterate through all non-trivial trees
    for k in 1:nt
        # Position of the first edge in the tree
        first = tree_edge_indices[k]

        # Total number of edges in the tree
        ne_tree = tree_edge_indices[k + 1] - first

        # Check if we have more than one edge in the tree (non-trivial tree)
        if ne_tree > 1
            # Determine if the tree is a star
            if is_star[k]
                # It is a non-trivial star and only the color of the hub is needed
                (leaf, hub) = reverse_bfs_orders[first]
                if hub < leaf
                    column_color_used[color[hub]] = true
                else
                    row_color_used[color[hub]] = true
                end
            else
                # It is not a star and both colors are needed during the decompression
                (i, j) = reverse_bfs_orders[first]
                v_col = min(i, j)
                v_row = max(i, j)
                row_color_used[color[v_row]] = true
                column_color_used[color[v_col]] = true
            end
        else
            nb_trivial_trees += 1
        end
    end

    # Process the trivial trees (if any)
    if nb_trivial_trees > 0
        # Color pairs (color of v_row, color of v_col) of the trivial trees whose root is not
        # determined yet: neither color is already required on its side.
        undecided_pairs = Tuple{Int,Int}[]

        for k in 1:nt
            first = tree_edge_indices[k]
            ne_tree = tree_edge_indices[k + 1] - first
            if ne_tree == 1
                (i, j) = reverse_bfs_orders[first]
                v_col = min(i, j)
                v_row = max(i, j)
                if column_color_used[color[v_col]]
                    # The vertex v_col is already an internal node in a non-trivial tree
                    reverse_bfs_orders[first] = (v_row, v_col)
                elseif row_color_used[color[v_row]]
                    # The vertex v_row is already an internal node in a non-trivial tree
                    reverse_bfs_orders[first] = (v_col, v_row)
                else
                    push!(undecided_pairs, (color[v_row], color[v_col]))
                end
            end
        end

        if !isempty(undecided_pairs)
            # Decide which colors to keep for the remaining trivial trees
            cover_trivial_structures!(
                row_color_used,
                column_color_used,
                undecided_pairs,
                bicoloring,
                postprocessing_minimizes,
            )

            # Every trivial tree now has an endpoint with a used color: make it the root
            for k in 1:nt
                first = tree_edge_indices[k]
                ne_tree = tree_edge_indices[k + 1] - first
                if ne_tree == 1
                    (i, j) = reverse_bfs_orders[first]
                    v_col = min(i, j)
                    v_row = max(i, j)
                    if column_color_used[color[v_col]]
                        reverse_bfs_orders[first] = (v_row, v_col)
                    else
                        reverse_bfs_orders[first] = (v_col, v_row)
                    end
                end
            end
        end
    end
    return nothing
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

## Choice of the colors kept for the trivial two-colored structures

#=
After the forced steps of the post-processing, each remaining trivial star or tree {i, j},
with i > j, still needs one of its two colors: the color of i on its side (the row side in a
bicoloring) or the color of j on its side (the column side). Choosing which colors to keep is a
vertex cover problem on a small "color graph", whose vertices are the colors (one copy per side
in a bicoloring) and whose edges are the distinct pairs of colors of the undecided structures.
Keeping a color costs one color whatever the number of structures it settles.

- In a bicoloring, the color graph is bipartite (row colors versus column colors), so a minimum
  vertex cover is computed exactly with Kőnig's theorem from a maximum matching. With
  `:row_colors` (resp. `:column_colors`), all column (resp. row) colors are kept instead, which
  never adds a color on the side to minimize.
- For a symmetric coloring, the color graph is general and the cover is chosen greedily, by
  repeatedly keeping the color that settles the largest number of remaining pairs.
=#

function cover_trivial_structures!(
    row_color_used::Vector{Bool},
    column_color_used::Vector{Bool},
    undecided_pairs::Vector{Tuple{Int,Int}},
    bicoloring::Bool,
    postprocessing_minimizes::Symbol,
    matching::AbstractMatchingAlgorithm=KuhnMatching(),
)
    if bicoloring
        if postprocessing_minimizes == :all_colors
            # The column colors are the "left" side of the bipartite color graph: on a single
            # undecided structure, the cover then keeps the column color, as the rule used
            # before did when nothing else decided.
            nb_colors = length(row_color_used)
            ptr, adj = _bipartite_adjacency(
                [(cj, ci) for (ci, cj) in undecided_pairs], nb_colors, nb_colors
            )
            bipartite_vertex_cover!(column_color_used, row_color_used, ptr, adj, matching)
        elseif postprocessing_minimizes == :row_colors
            for (_, cj) in undecided_pairs
                column_color_used[cj] = true
            end
        elseif postprocessing_minimizes == :column_colors
            for (ci, _) in undecided_pairs
                row_color_used[ci] = true
            end
        else
            error(
                "The value postprocessing_minimizes = :$postprocessing_minimizes is not supported.",
            )
        end
    else
        # row_color_used and column_color_used point to the same memory
        greedy_cover!(row_color_used, undecided_pairs)
    end
    return nothing
end

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

function greedy_cover!(color_used::Vector{Bool}, pairs::Vector{Tuple{Int,Int}})
    nb_colors = length(color_used)
    # undirected simple graph on the colors, each pair stored in both directions
    sym_pairs = Tuple{Int,Int}[]
    for (a, b) in pairs
        push!(sym_pairs, (a, b))
        push!(sym_pairs, (b, a))
    end
    sym_pairs = unique!(sort!(sym_pairs))
    ptr, adj = _pair_adjacency(sym_pairs, nb_colors)

    # dynamic degree = number of uncovered pairs incident to the color
    degree = [ptr[c + 1] - ptr[c] for c in 1:nb_colors]
    dmax = maximum(degree; init=0)
    # buckets with lazy deletion: a color may appear in several buckets, only the one
    # matching its current degree is valid
    buckets = [Int[] for _ in 1:dmax]
    for c in nb_colors:-1:1
        degree[c] > 0 && push!(buckets[degree[c]], c)
    end
    d = dmax
    while d > 0
        if isempty(buckets[d])
            d -= 1
            continue
        end
        c = pop!(buckets[d])
        (color_used[c] || degree[c] != d) && continue
        # keep the color c, which covers all its remaining pairs
        color_used[c] = true
        degree[c] = 0
        for p in ptr[c]:(ptr[c + 1] - 1)
            c2 = adj[p]
            if !color_used[c2]
                degree[c2] -= 1
                degree[c2] > 0 && push!(buckets[degree[c2]], c2)
            end
        end
    end
    return nothing
end
