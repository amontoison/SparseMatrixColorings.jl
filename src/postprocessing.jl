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
