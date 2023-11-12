## Work on sizing
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
    ' ' => 9,
    '*' => 9, 
)
clustal_colormap = tuple.([:blue, :red, :magenta, :green, :pink, :orange, :yellow, :cyan, :white], 0.5)

msa_pf62 = MIToS.MSA.read("pf00062.stockholm.gz",Stockholm)
window = 50
seqs = 30
font_size = 12
cell_size = 14
padding = 200

msa1 = permutedims(Char.(msa_pf62[1:seqs,:]))

mul1(x,s) = (x-1)*s +1

function msaplot1!(ax, msa_matrix, len, sl; color_map = clustal_colormap, res_levels = clustal_levels, scale = 1)
    mat = @lift(msa_matrix[mul1($(sl.value),scale):min(mul1($(sl.value)+len,scale), end),:])
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
function positions!(ax, len, max_pos, sl; scale = 1)
    xs = @lift(string.(mul1($(sl.value),scale):mul1(min($(sl.value)+len,scale), max_pos)))
    text!(ax, 1:mul1(len+1,scale), repeat([1],mul1(len+1,scale)), text=xs,
          align = (:center, :center),
          rotation=-pi/2,
          )
    hidedecorations!(ax)
    xlims!(0.5,mul1(len+1, scale) + .5)
end

f = Figure(fontsize = font_size, resolution = (3*window*cell_size + padding, (3*seqs)*cell_size))
sl = Makie.Slider(f[3, 1:3], range = 1:size(msa1,1)-window, startvalue = 1)
msaplot1!(Axis(f[4,1]; yreversed=true,), msa1, window, sl)
positions!(Axis(f[5,1],height = 50), window, size(msa1,1),sl)

## dna
using FASTX
using BioSequences
dna1 = first(FASTAReader(open("s1.fasta")))
dna = sequence(LongDNA{4},dna1)
t1 = translate(dna)
t2 = translate(dna[2:end-2])
t3 = translate(dna[3:end-1])

function v2char(x, sep="") collect.(join(string.(x), sep)) end

dna_mat = [v2char(dna) [v2char(t1,"  ") ;[' ', ' ']] [[' '] ;v2char(t2,"  ") ;[ ' ',' ',' ',' ']] [[' ',' '] ;v2char(t3,"  ") ;[ ' ',' ',' ']]]
msaplot1!(Axis(f[1,1:3]; yreversed=true,height=80), dna_mat, window, sl; scale = 3)
positions!(Axis(f[2,1:3],height = 50), window, size(msa1,1),sl; scale = 3)
