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
    heatmap!(ax, l1, colormap = color_map, colorrange=(1,length(color_map)), inspectable = false)
    t1 = @lift(vec(string.($mat)))
    p1 = @lift(vec(Point2f.(Tuple.(CartesianIndices($mat)))))
    text!(ax, t1,
          position = p1,
          align = (:center, :center),
          )
    hidedecorations!(ax)
end

function o_plot_mat!(ax, mat)
    t1 = @lift(vec($mat))
    p1 = @lift(vec(Point2f.(Tuple.(CartesianIndices($mat)))))
    ## text!(ax, t1, position = p1, align=(:center, :center), rotation=-pi/2)
    text!(ax, p1; text = t1, align=(:center, :center), rotation=-pi/2)
    xlims!(0.5, size(mat[])[1]+.5)
    hidexdecorations!(ax)
end

# function o_plot_dna1!(ax, mat, dna_pos; color_map = clustal_colormap, res_levels = clustal_levels)
#     frame1_true = @lift($dna_pos[:,end] .== "1")
#     frame1 = @lift(parse.(Int,$dna_pos[$frame1_true,1]) .- parse(Int,$dna_pos[1,1]) .+ .5)
#     l1 = @lift([res_levels[m] for m in $mat])
#     vlines!(ax, frame1, color = (:black,1.))
#     heatmap!(ax, l1, colormap = color_map, colorrange=(1,length(color_map)))
#     t1 = @lift(vec(string.($mat)))
#     p1 = @lift(vec(Point2f.(Tuple.(CartesianIndices($mat)))))
#     text!(ax, t1,
#           position = p1,
#           align = (:center, :center), # inspector_label
#           inspector_label = (self, i, p) ->  "HEP!", ## Does not work for test and volume  https://docs.makie.org/stable/explanations/inspector/index.html
#           )
#     hidedecorations!(ax)
# end


function o_plot_dna!(ax, mat, dna_pos; color_map = clustal_colormap, res_levels = clustal_levels)
    frame1_true = @lift($dna_pos[:,end] .== "1")
    frame1 = @lift(parse.(Int,$dna_pos[$frame1_true,1]) .- parse(Int,$dna_pos[1,1]) .+ .5)
    l1 = @lift([res_levels[m] for m in $mat])
    vlines!(ax, frame1, color = (:black,2))
    heatmap!(ax, l1, colormap = color_map, colorrange=(1,length(color_map)), inspectable = false)
    t1 = @lift((vec($mat))) # uppercase.
    p1 = @lift(vec(Point2f.(Tuple.(CartesianIndices($mat)))))
    scatter!(ax, p1, marker = t1,
             align = (:center, :center),
             color = :black,
             inspector_label = (self, i, p) -> i <= size(mat[],1) ? "$(uppercase.(mat[][i,1]))$(dna_pos[][i])" : ""  , ## Does not work for test and volume  https://docs.makie.org/stable/explanations/inspector/index.html  "i: $i" ##
          )
    hidedecorations!(ax)
end

function o_plot_dna2!(ax, mat, dna_pos; color_map = clustal_colormap, res_levels = clustal_levels)
    ## color by exon structure: col 2 of dna_pos
    ## TODO: include coloring of reading-frame?
    frame1_true = @lift($dna_pos[:,end] .== "1")
    frame1 = @lift(parse.(Int,$dna_pos[$frame1_true,1]) .- parse(Int,$dna_pos[1,1]) .+ .5)
    l1 = @lift([res_levels[m] for m in $dna_pos[:,2]]) ## intron, exon
    l2 = @lift(repeat([2], size($dna_pos,1)))
    l4 = @lift(hcat($l1, $l2, $l2, $l2))
    vlines!(ax, frame1, color = (:black,2))
    heatmap!(ax, l4, colormap = color_map, colorrange=(1,length(color_map)), inspectable = false)
    t1 = @lift((vec($mat))) # uppercase.
    p1 = @lift(vec(Point2f.(Tuple.(CartesianIndices($mat)))))
    scatter!(ax, p1, marker = t1,
             align = (:center, :center),
             color = :black,
             inspector_label = (self, i, p) -> i <= size(mat[],1) ? "$(uppercase.(mat[][i,1]))$(dna_pos[][i])" : ""  , ## Does not work for test and volume  https://docs.makie.org/stable/explanations/inspector/index.html  "i: $i" ##
          )
    hidedecorations!(ax)
end

## TODO: module, intron plot:
## show modules (of homologues) and mark introns by vertical lines
## This plot follows MSA plot

## TODO click sequences out (and back in of the MSA (pop to separate track)
