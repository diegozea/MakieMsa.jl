function positions!(ax, xs)
    npos = length(xs[])
    text!(ax, 1:npos, repeat([1],npos), text=xs,
          align = (:center, :center),
          rotation=-pi/2,
          )
    hidedecorations!(ax)
    xlims!(0.5,npos+.5)
end

function dna_positions!(ax, dna_pos, aa_pos)
    npos = length(dna_pos[])
    text!(ax, 1:npos, repeat([1],npos), text=dna_pos,
          align = (:center, :center),
          rotation=-pi/2,
          )
    hidedecorations!(ax)
    xlims!(0.5,npos+.5)
end


function msaplot!(ax, mat;  color_map = clustal_colormap, res_levels = clustal_levels)
    l1 = @lift([res_levels[m] for m in $mat])
    heatmap!(ax, l1, colormap = color_map, colorrange=(1,length(color_map)))
    t1 = @lift(vec(string.($mat)))
    p1 = @lift(vec(Point2f.(Tuple.(CartesianIndices($mat)))))
    text!(ax, t1,
          position = p1,
          align = (:center, :center),
          )
    hidedecorations!(ax)
end

function plot_mat!(ax, mat)
    t1 = @lift(vec($mat))
    p1 = @lift(vec(Point2f.(Tuple.(CartesianIndices($mat)))))    
    ## text!(ax, t1, position = p1, align=(:center, :center), rotation=-pi/2)
    text!(ax, p1; text = t1, align=(:center, :center), rotation=-pi/2)
    xlims!(0.5, size(mat[])[1]+.5)
    hidexdecorations!(ax)
end