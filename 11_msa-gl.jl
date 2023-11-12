using WGLMakie
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
sl = Makie.Slider(f[2, 1], range = 1:5, startvalue = 1)

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
msaplot5!(Axis(f[1,1]; yreversed=true,), msa_matrix , colormap,sl)

## Now add position!
function positions!(ax, len, max_pos, sl)
    xs = @lift(string.($(sl.value):min($(sl.value)+len, max_pos)))
    text!(ax, 1:len+1, repeat([1],len+1), text=xs,
          align = (:center, :center),          
          )    
end

positions!(Axis(f[3,1]), 6,20,sl)

function msaplot6!(ax, msa_matrix, colors, len, sl)
    mat = @lift(msa_matrix[$(sl.value):min($(sl.value)+len, end),:])
    l1 = @lift([levels[m] for m in $mat])
    heatmap!(ax, l1, colormap = colors, colorrange=(1,5))
    t1 = @lift(vec(string.($mat)))
    p1 = @lift(vec(Point2f.(Tuple.(CartesianIndices($mat)))))
    text!(ax, t1,
          position = p1,
          align = (:center, :center),          
          )    
end

## Add bigger one + positions
msaplot6!(Axis(f[4,1:2]; yreversed=true,), msa_matrix , colormap, 12, sl)
positions!(Axis(f[5,1:2]), 12,20,sl)
