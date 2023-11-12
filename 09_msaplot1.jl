# https://discourse.julialang.org/t/makie-multi-sequence-alignment-scatterplot/77932
using GLMakie

msa_matrix = permutedims([
    '-' 'G' 'G' '-' 'C' 'T' 'T' 'G' 'C' 'T' 'T' 'A' 'T' 'T' 'G' 'T' '-' '-' 'G' 'T'
    '-' 'G' 'G' 'C' 'T' 'T' 'G' 'C' 'T' 'T' 'A' 'T' 'T' 'G' 'T' 'G' 'T' 'G' 'G' 'T'
    'C' 'G' 'G' '-' 'T' 'T' 'G' 'C' 'T' 'T' 'A' 'T' 'T' 'G' 'T' 'G' '-' '-' '-' 'T'
    'C' 'G' 'G' '-' 'T' 'T' 'G' 'C' 'T' 'T' 'A' 'T' 'T' 'G' 'T' 'G' '-' '-' 'T' '-'
])

levels = Dict(
    'A' => 1,
    'C' => 2,
    'G' => 3,
    'T' => 4,
    '-' => 5,
)
colormap = tuple.([:yellow, :blue, :green, :red, :magenta], 0.5)

f = Figure()
ax = Axis(f[1,1])
hidespines!(ax)

heatmap!(ax, [levels[m] for m in msa_matrix], colormap = colormap)
text!(ax, vec(string.(msa_matrix)),
    position = vec(Point2f.(Tuple.(CartesianIndices(msa_matrix)))),
    align = (:center, :center),

)

ax2 = Axis(f[2,1],     xgridvisible = false,
    ygridvisible = false,
)


function msaplot1!(fig_pos, mat, colors)
    ax = Axis(fig_pos,
              xgridvisible = false,
              ygridvisible = false,
              xticks = 1:size(mat, 1),
              yticks = (1:size(mat, 2), ["a", "b", "c", "d"]),
              )
    heatmap!(ax, [levels[m] for m in mat], colormap = colors)
    text!(ax, vec(string.(mat)),
          position = vec(Point2f.(Tuple.(CartesianIndices(mat)))),
          align = (:center, :center),          
          )    
end

msaplot1!(f[2,1], msa_matrix, colormap)
msaplot1!(f[3,1], msa_matrix, colormap)
msaplot1!(f[4,1], msa_matrix, colormap)

function msaplot2!(fig_pos, mat, colors, start, len)
    mat = mat[start:min(start+len, end),:]
    ax = Axis(fig_pos,
              xgridvisible = false,
              ygridvisible = false,
              xticks = start:start+len,
              yticks = (1:size(mat, 2), ["a", "b", "c", "d"]),
              )
    empty!(ax)
    heatmap!(ax, [levels[m] for m in mat], colormap = colors)
    text!(ax, vec(string.(mat)),
          position = vec(Point2f.(Tuple.(CartesianIndices(mat)))),
          align = (:center, :center),          
          )    
end

msaplot2!(f[5,1], msa_matrix, colormap, 5,5)
msaplot2!(f[6,1], msa_matrix, colormap, 7,5)

sl = Slider(f[7, 1], range = 1:5, startvalue = 1)



function msaplot3!(fig_pos, mat, colors, slider, len)
    msaplot2!(fig_pos, mat, colors, slider.val, len)
end

## Restart
f = Figure()

function msaplot4!(ax, mat, colors)
    l1 = @lift([levels[m] for m in $mat])
    heatmap!(ax, l1, colormap = colors, colorrange=(1,5))
    t1 = @lift(vec(string.($mat)))
    p1 = @lift(vec(Point2f.(Tuple.(CartesianIndices($mat)))))
    text!(ax, t1,
          position = p1,
          align = (:center, :center),          
          )    
end

ax1 = Axis(f[1,1])
sl = Slider(f[2, 1], range = 1:5, startvalue = 1)
mat = @lift(mat = msa_matrix[$(sl.value):min($(sl.value)+6, end),:])
            
msaplot4!(ax1, mat, colormap)

function msaplot5!(ax, msa_matrix, colors, sl)
    mat = @lift(msa_matrix[$(sl.value):min($(sl.value)+6, end),:])
    l1 = @lift([levels[m] for m in $mat])
    heatmap!(ax, l1, colormap = colors, colorrange=(1,5))
    t1 = @lift(vec(string.($mat)))
    p1 = @lift(vec(Point2f.(Tuple.(CartesianIndices($mat)))))
    text!(ax, t1,
          position = p1,
          align = (:center, :center),          
          )    
end

msaplot5!(Axis(f[3,1]), msa_matrix , colormap,sl)

function msaplot6!(fig_pos, msa_matrix, colors, sl)
    x_ticks = sl.value.val:sl.value.val+6
    #center_val = string(sl.value.val + 3)
    ax = Axis(fig_pos,
              xgridvisible = false,
              ygridvisible = false,
              xticks = x_ticks,
              yticks = (1:4, ["a", "b", "c", "d"]),
              # title=center_val,
              )
    
    mat = @lift(msa_matrix[$(sl.value):min($(sl.value)+6, end),:])
    l1 = @lift([levels[m] for m in $mat])
    heatmap!(ax, l1, colormap = colors, colorrange=(1,5))
    t1 = @lift(vec(string.($mat)))
    p1 = @lift(vec(Point2f.(Tuple.(CartesianIndices($mat)))))
    text!(ax, t1,
          position = p1,
          align = (:center, :center),          
          )    
end

msaplot6!(Axis(f[4,1]), msa_matrix , colormap, sl)
