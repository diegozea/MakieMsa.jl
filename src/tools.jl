function v2char(x, sep="") collect.(join(string.(x), sep)) end
function div1(x,y) div(x-1,y)+1 end


## Convert DNA pos to codon pos:
function dna2codon(dna_start)
    div1(dna_start,3)
end
