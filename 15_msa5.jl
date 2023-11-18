## Work on sizing
using GLMakie
using WGLMakie
using DataStructures
using MIToS
using MIToS.MSA, MIToS.Pfam
## dna
using FASTX
using BioSequences


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


function positions1!(ax, len, max_pos, sl)
    xs = @lift(string.($(sl.value):min($(sl.value)+len, max_pos)))
    text!(ax, 1:len+1, repeat([1],len+1), text=xs,
          align = (:center, :center),
          rotation=-pi/2,
          )
    hidedecorations!(ax)
    xlims!(0.5,len+1.5)
end

function v2char(x, sep="") collect.(join(string.(x), sep)) end

function div1(x,y) div(x-1,y)+1 end

## Convert DNA pos to msa pos:
function dna2msa(dna_start)
    ## demo: just div 3
    max(div(dna_start,3),1)
end
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

const msa_pf62 = MIToS.MSA.read("pf00062.stockholm.gz",Stockholm)
const dna1 = first(FASTAReader(open("s1.fasta")))
const dna = sequence(LongDNA{4},dna1)
const t1 = translate(dna)
const t2 = translate(dna[2:end-2])
const t3 = translate(dna[3:end-1])
const dna_mat = [v2char(dna) [v2char(t1,"  ") ;[' ', ' ']] [[' '] ;v2char(t2,"  ") ;[ ' ',' ',' ',' ']] [[' ',' '] ;v2char(t3,"  ") ;[ ' ',' ',' ']]]



const msa_window = 50
const dna_window = 3*msa_window
const msa_seqs = 20
const font_size = 12
const cell_size = 14
const padding = 200

msa1 = permutedims(Char.(msa_pf62[1:msa_seqs,:]))

dna_length = length(dna)
msa_length = size(msa1,1)

f = Figure(fontsize = font_size, resolution = (dna_window*cell_size + padding, 3*msa_seqs*cell_size))
sl = Makie.Slider(f[3, 1:3], range = 1:length(dna), startvalue = 1) ## Slider

## Compute positions
dna_pos = @lift(string.($(sl.value):min($(sl.value)+dna_window, dna_length)))
msa_pos = @lift(string.(dna2msa($(sl.value)):min(dna2msa($(sl.value)) + msa_window, msa_length)))
aa_pos = @lift(div1.(parse.(Int, $dna_pos),3)[1:length($dna_pos)])

positions!(Axis(f[1,1:3],height = 50), dna_pos) ## DNA position
positions!(Axis(f[5,1:2],height = 50), msa_pos) ## MSA position

## Compute matrices to show
dna_show = @lift(dna_mat[parse.(Int,$dna_pos),:])
msa_show = @lift(msa1[parse.(Int,$msa_pos),:])

msaplot!(Axis(f[2,1:3]; yreversed=true,height=80), dna_show) ## DNA
msaplot!(Axis(f[4,1:2]; yreversed=true), msa_show) ## MSA

## GLMakie.activate!()
## WGLMakie.activate!()

## https://diegozea.github.io/MIToS.jl/latest/MSA/#Column-and-sequence-mappings
## a2 = read("https://raw.githubusercontent.com/diegozea/MIToS.jl/master/test/data/PF09645_full.stockholm", Stockholm, generatemapping=true, useidcoordinates=false, deletefullgaps=false)
## getsequencemapping(a2, "C3N734_SULIY/1-95")
## getsequence(a2,4)
## stringsequence(a2,1) 
