## Add clustal colors
using GLMakie
using DataStructures
using MIToS
using MIToS.MSA, MIToS.Pfam


## Clustal colors: https://www.jalview.org/help/html/colourSchemes/clustal.html
clustal_levels = OrderedDict(
    'A' => 1,
    'B' => 9,
    'C' => 6,
    'D' => 3,
    'E' => 3,
    'F' => 1,
    'G' => 6,
    'H' => 8,
    'I' => 1,
    'J' => 9,
    'K' => 2,
    'L' => 1,
    'M' => 1,
    'N' => 4,
    'O' => 9,
    'P' => 7,
    'Q' => 4,
    'R' => 2,
    'S' => 4,
    'T' => 4,
    'U' => 9,
    'V' => 1,
    'W' => 1,
    'X' => 9,
    'Y' => 8,
    'Z' => 9,
    '-' => 9,
    '.' => 9,
)
clustal_colormap = tuple.([:blue, :red, :magenta, :green, :pink, :orange, :yellow, :cyan, :white], 0.5)

msa_pf62 = MIToS.MSA.read("pf00062.stockholm.gz",Stockholm)
msa1 = permutedims(Char.(msa_pf62[1:30,:]))
window = 50

function msaplot1!(ax, msa_matrix, len, sl; color_map = clustal_colormap, res_levels = clustal_levels)
    mat = @lift(msa_matrix[$(sl.value):min($(sl.value)+len, end),:])
    l1 = @lift([res_levels[m] for m in $mat])
    heatmap!(ax, l1, colormap = color_map, colorrange=(1,length(color_map)))
    t1 = @lift(vec(string.($mat)))
    p1 = @lift(vec(Point2f.(Tuple.(CartesianIndices($mat)))))
    text!(ax, t1,
          position = p1,
          align = (:center, :center),          
          )    
end
function positions!(ax, len, max_pos, sl)
    xs = @lift(string.($(sl.value):min($(sl.value)+len, max_pos)))
    text!(ax, 1:len+1, repeat([1],len+1), text=xs,
          align = (:center, :center),
          rotation=-pi/2,
          )    
end

f = Figure()
sl = Makie.Slider(f[2, 1:2], range = 1:100, startvalue = 1)
msaplot1!(Axis(f[3,1:2]; yreversed=true,), msa1, 30, sl)
positions!(Axis(f[4,1:2]), window, size(msa1,1),sl)
