"""
    dna_location
    Compute coordinates on DNA sequence given location strin indication exon structure
    Input:
    dna: a DNA sequence
    location: a location string. Eg join(1..10,15..20)
"""
function dna_positions(dna, location)
    exons = exons_from_location(location)
    
end


Base.range(m::RegexMatch) = m.offset .+ (0:length(m.match)-1) ## type pirate! From https://discourse.julialang.org/t/find-index-of-all-occurences-of-a-string/23044/8

substring(s::String,v::UnitRange) = s[v]
function substring(s::String,r::Regex)
    idx = range.(eachmatch(r, s))
    substring.(s, idx)
end

function exons_from_location(location)
    reg = r"(\d+\D+\d+)" ## 123..234,324..456
    exon_strings = substring(location, reg) ## ["123..345", "324..456"]
    map(x-> parse.(Int, x), split.(exon_strings, r"\.+")) ## [[123,234],[324,456]]
end
    
