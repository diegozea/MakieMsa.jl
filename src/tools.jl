function v2char(x, sep=""; tolower=false)
    x1 = join(string.(x), sep)
    if tolower
        x1 = lowercase(x1)
    end
    collect.(x1)
end
function div1(x,y) div(x-1,y)+1 end


## Convert DNA pos to codon pos:
function dna2codon(dna_start)
    div1(dna_start,3)
end
