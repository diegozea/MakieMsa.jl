## Work on sizing
#using WGLMakie
using GLMakie
using DataStructures
using MIToS
using MIToS.MSA, MIToS.Pfam
## dna
using FASTX
using BioSequences

include("src/consts.jl")
include("src/tools.jl")
include("src/plots.jl")
include("src/dna_location.jl")
using DataFrames

# ## Convert DNA pos to msa pos:
# function dna2msa(dna_start)
#     ## demo: just div 3
#     max(div(dna_start,3),1)
# end


const msa_pf62 = MIToS.MSA.read("pf00062.stockholm.gz",Stockholm)
dna1 = first(FASTAReader(open("s1.fasta")))
dna2 = sequence(LongDNA{4},dna1)
dna3 = string(dna2)
intron1 = join(rand(['A','C','G','T'],100))
intron2 = join(rand(['A','C','G','T'],100))
intron3 = join(rand(['A','C','G','T'],49))
gene_model = "1..50,151..200,301..551"
dna4 = dna3[1:50] * intron1 * dna3[51:100] * intron2 * dna3[101:end] * intron3
dna = LongDNA{4}(dna4)

dna_pos = dna_positions(gene_model, [1,length(dna4)])

t1 = translate(dna)
t2 = translate(dna[2:end-2])
t3 = translate(dna[3:end-1])
dna_mat = [v2char(dna; tolower=true) [v2char(t1,"  ") ;[' ', ' ']] [[' '] ;v2char(t2,"  ") ;[ ' ',' ',' ',' ']] [[' ',' '] ;v2char(t3,"  ") ;[ ' ',' ',' ']]]



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
dna_pos_mat = replace(string.(Matrix(dna_pos)), "missing" => "", "exon" => "e", "intron" => "i", "UTR" => "u")
o_dna_range = @lift($(sl.value):min($(sl.value)+dna_window, dna_length))
o_dna_pos = @lift(dna_pos_mat[$o_dna_range,:])
## msa_pos = @lift(string.(dna2codon($(sl.value)):min(dna2codon($(sl.value)) + msa_window, msa_length)))
o_codon_pos = @lift(minimum(skipmissing(dna_pos.codon_number[$o_dna_range])))
msa_pos = @lift(string.($o_codon_pos: min($o_codon_pos + msa_window, msa_length)))

o_plot_mat!(Axis(f[1,1:3], height = 50*4+10, yticks = (1:8, names(dna_pos))), o_dna_pos) ## DNA position
positions!(Axis(f[5,1:2],height = 50), msa_pos) ## MSA position

## Compute matrices to show
dna_show = @lift(dna_mat[$o_dna_range,:])
msa_show = @lift(msa1[parse.(Int,$msa_pos),:])

o_plot_dna!(Axis(f[2,1:3]; yreversed=true,height=80), dna_show, o_dna_pos) ## DNA
msaplot!(Axis(f[4,1:2]; yreversed=true), msa_show) ## MSA

## DataInspector(Axis(f[2,1:3]))

GLMakie.activate!()
## WGLMakie.activate!()

## https://diegozea.github.io/MIToS.jl/latest/MSA/#Column-and-sequence-mappings
## a2 = read("https://raw.githubusercontent.com/diegozea/MIToS.jl/master/test/data/PF09645_full.stockholm", Stockholm, generatemapping=true, useidcoordinates=false, deletefullgaps=false)
## getsequencemapping(a2, "C3N734_SULIY/1-95")
## getsequence(a2,4)
## stringsequence(a2,1) 

## Serving from Jupyter:
## https://github.com/SimonDanisch/JSServe.jl
### example: https://github.com/SimonDanisch/SmartHomy/blob/master/web_app.jl
