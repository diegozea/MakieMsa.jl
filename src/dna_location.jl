"""
    dna_positions
    Compute DataFrame with coordinates in genomic DNA and exon given exon location string (join) and linear range of DNA postions.
    exon locations must be completely inside the dna range.
    Input:
    exon_location. Location string. Eg "1..200" or join(10..20,30..40)
    dna_range: range of the dna positions. Eg [1,300]
    Output:
    DataFrame with dna_pos  state   exon_number  exon_pos  intron_number  intron_pos  codon_number  codon_pos
    dna_pos: position in genomic dna
    state: exon, intron or UTR (before and after all exons)
    exon_number: counts the exons from 1. Missing outside of exons
    exon_pos: position in spliced exon increasing from 1. Missing outside of exons
    intron_number: counts introns from 1. Missing outside of introns
    intron_pos: position in current intron. Reset to 1 in each intron.
    codon_number: codon number in exon. This will correspond to position in protein sequence. Missing outside of exons.
    codon_pos: position in codon. Counts 1, 2, 3 in each codon. Missing outside of exons. Must end as 3 (end in frame)
"""
function dna_positions(exon_location, dna_range)
    @assert( minimum(dna_range) > 0 && dna_range[2] > dna_range[1] )
    exons = exons_from_location(exon_location)
    @assert length(exons) >= 1
    dna_pos = Int[]
    state = String[]
    exon_number = Union{Int, Missing}[]
    exon_pos = Union{Int, Missing}[]
    intron_number = Union{Int, Missing}[]
    intron_pos = Union{Int, Missing}[] ## reset for each intron
    codon_number = Union{Int, Missing}[]
    codon_pos = Union{Int, Missing}[] ## reset for each codon, ie local readingframe
    local d_pos = dna_range[1]
    local current_state = "UTR"
    current_exon = popfirst!(exons)
    local e_num = 0
    local e_pos = 1
    local i_num = 0
    local i_pos = 1
    local c_num = 1
    local c_pos = 1
    while d_pos <= dna_range[2]
        push!(dna_pos, d_pos)
        if d_pos == current_exon[1]  ## enter exon
            current_state = "exon"
            e_num += 1 ## next exon
            i_pos = 1 ## reset intron pos
        end
        push!(state, current_state)
        if current_state == "exon"
            push!(exon_pos, e_pos)
            e_pos += 1
            push!(exon_number, e_num)
        else
            push!(exon_pos, missing)
            push!(exon_number, missing)
        end
        if current_state == "intron"
            push!(intron_pos, i_pos)
            i_pos += 1
            push!(intron_number, i_num)
        else
            push!(intron_pos, missing)
            push!(intron_number, missing)
        end
        d_pos += 1
        if d_pos == current_exon[2]+1
            if length(exons) > 0
                current_state = "intron"
                current_exon = popfirst!(exons)
                i_num += 1
            else
                current_state = "UTR"
                current_exon = dna_range[2] .+ dna_range ## will never reach
            end
        end
    end
    codon_number = div1.(exon_pos,3)
    codon_pos= mod1.(exon_pos,3)
    @assert collect(skipmissing(codon_pos))[end] == 3 ## ends in frame
    DataFrame(;dna_pos, state, exon_number, exon_pos, intron_number, intron_pos, codon_number, codon_pos)
end


Base.range(m::RegexMatch) = m.offset .+ (0:length(m.match)-1) ## type pirate! From https://discourse.julialang.org/t/find-index-of-all-occurences-of-a-string/23044/8

substring(s::String,v::UnitRange) = s[v]
function substring(s::String,r::Regex)
    idx = range.(eachmatch(r, s))
    substring.(s, idx)
end

"""
    parse_location
    Parse a location string to a vector of 2-vectors or the endpoints.
    For now, only fwd strand is supported and only actual enpoints (not 1..>123)
    Examples:
    parse_location("1..200") == [[1,200]]
    parse_location(join("1..200,300..400")) == [[1,200], [300,400]]
"""
function parse_location(location)
    reg = r"(\d+\D+\d+)" ## 123..234,324..456
    exon_strings = substring(location, reg) ## ["123..345", "324..456"]
    map(x-> parse.(Int, x), split.(exon_strings, r"\.+")) ## [[123,234],[324,456]]
end
    
function div1(x,y) div(x-1,y)+1 end
import Base.mod1 ## type pirate!
mod1(::Missing,x) = missing
