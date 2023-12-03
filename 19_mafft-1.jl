## Base on mafft
using GLMakie
using FASTX
using BioSequences
using DataFrames
using GenomicAnnotations ## for reading genbankfiles 
include("src/mafft.jl")

## align
msa_file = mafft("data/aa_1.fas")

msa1 = read_msa(msa_file)

mat1 = msa_matrix(msa1, reverse_order=false)

f = Figure() ## fontsize = font_size, resolution = (dna_window*cell_size + padding, 3*msa_seqs*cell_size))

plot_charmatrix!(Axis(f[1,1]), mat1)

## reading DNA:
## using GenomicAnnotations
d1 = readgbk("data/NZGN_EFN1HTZPWW_lflank_500_start_91855_stop_92956_rflank_500.gb")
g1 = @genes(d1, CDS)
## s1 = GenomicAnnotations.sequence(d1[1].genes[1])
s1 = GenomicAnnotations.sequence(g1[1])
l1 = locus(g1[1])
## from https://github.com/BioJulia/GenomicAnnotations.jl/blob/573e0226d5473cec7743fbed837124346bbd17ff/src/record.jl
# struct Locus
#     position::UnitRange{Int}
#     strand::Char ## '+' or '-'
#     complete_left::Bool
#     complete_right::Bool
#     order::Vector{UnitRange{Int}}
#     join::Bool
# end
# function sequence(gene::AbstractGene; translate = false, preserve_alternate_start = false)
#     if locus(gene).strand == '-'
#         s = reverse_complement(parent(gene).sequence[locus(gene).position])
#     else
#         s = parent(gene).sequence[locus(gene).position]
#     end

## TODO: Handle locus on reverse strand: convert before editing, and convert after editing
function rev_range(range, genome_length)
    t1 = reverse(genome_length .- extrema(range) .+1)
    t1[1]:t1[2]
end
function reverse_locus(l1, d1)
    if l1.strand == '+'
        @info "locus already forward: $(l1)"
        return (l1, d1)
    end
    if !(l1.complete_left & l1.complete_right)
        error("Can not reverse partial location")
    end
    ## TODO: not join
    genome_length = length(d1.sequence)
    l2 = Locus(
        rev_range(l1.position, genome_length),
        '+',
        true,
        true,
        reverse(rev_range.(l1.order, genome_length)),
        l1.join,
    )
    (l2, BioSequences.reverse_complement(d1.sequence))
end

loc, dna1  = reverse_locus(l1, d1[1]) ## Always safe to call

spliced_dna = LongDNA{4}(join(map(x->dna1[x], loc.order)))
spliced_aa = BioSequences.translate(spliced_dna)

## from here basically 17_intron2.jl
using GLMakie
using DataStructures
using FASTX
using BioSequences

include("src/consts.jl")
include("src/tools.jl")
include("src/plots.jl")
include("src/dna_location.jl")
using DataFrames
gene_model = string(loc) ## "join(501..587,700..1085,1233..1602)"

dna = dna1[1:div(length(dna1) ,3)*3] ## truncate to codon
dna_pos = dna_positions(gene_model, [1,length(dna)])
t1 = translate(dna)
t2 = translate(dna[2:end-2])
t3 = translate(dna[3:end-1])
dna_mat = [v2char(dna; tolower=true) [v2char(t1,"  ") ;[' ', ' ']] [[' '] ;v2char(t2,"  ") ;[ ' ',' ',' ',' ']] [[' ',' '] ;v2char(t3,"  ") ;[ ' ',' ',' ']]]


dna_width::Int = 3
msa_width::Int = 3
msa_window::Int = 150
dna_window::Int = 1*msa_window
msa_seqs::Int = 20
font_size::Int = 12
cell_size::Int = 14
padding::Int = 200

msa1 = mat1  ## was permutedims(Char.(msa_pf62[1:msa_seqs,:])) now msa_matrix(read_msa(msaFile), reverse_order=true)

dna_length = length(dna)
msa_length = size(msa1,1)


f = Figure(fontsize = font_size, resolution = (dna_window*cell_size + padding, 3*msa_seqs*cell_size))
sl = Makie.Slider(f[3, 1:dna_width], range = 1:length(dna)-dna_window, startvalue = 1) ## Slider
## Compute positions
dna_pos_mat = replace(string.(Matrix(dna_pos)), "missing" => "", "exon" => "e", "intron" => "i", "UTR" => "u")
o_dna_range = @lift($(sl.value):min($(sl.value)+dna_window, dna_length))
o_dna_pos = @lift(dna_pos_mat[$o_dna_range,:])
## msa_pos = @lift(string.(dna2codon($(sl.value)):min(dna2codon($(sl.value)) + msa_window, msa_length)))
## o_codon_pos = @lift(minimum(coalesce(dna_pos.codon_number[$o_dna_range]..., msa_length)))
function first_after(v,p) ## first non-missing in v after pos p
    ## dna_pos.codon_number[findfirst(.!ismissing.(dna_pos.codon_number[331:end]))+331]
    v[findfirst(.!ismissing.(v[p:end]))+p]
end
o_codon_pos = @lift(first_after(dna_pos.codon_number, $(o_dna_range)[1]))
msa_pos = @lift(string.($o_codon_pos: min($o_codon_pos + msa_window, msa_length)))

o_plot_mat!(Axis(f[1,1:dna_width], height = 50*4+10, yticks = (1:8, names(dna_pos))), o_dna_pos) ## DNA position
positions!(Axis(f[5,1:msa_width],height = 50), msa_pos) ## MSA position

# Compute matrices to show 
dna_show = @lift(dna_mat[$o_dna_range,:])
msa_show = @lift(msa1[parse.(Int,$msa_pos),:])


## o_plot_dna!(Axis(f[2,1:dna_width]; yreversed=true,height=80), dna_show, o_dna_pos) ## DNA
o_plot_dna2!(Axis(f[2,1:dna_width]; yreversed=true,height=80), dna_show, o_dna_pos, color_map = exon_colormap, res_levels = exon_levels ) ## DNA
msaplot!(Axis(f[4,1:msa_width]; yreversed=true), msa_show) ## MSA

DataInspector(Axis(f[2,1:dna_width]))

GLMakie.activate!()
