using Test


@testset "dna_location" begin
    @test parse_location("join(123..234,324..456)") == [[123,234],[324,456]]
    @test parse_location("1..200") == [[1,200]]    
end

dna_pos = dna_positions("2..3,6..8,10..13", [1,15])
# julia> dna_positions("2..3,6..8,10..13", [1,15])
# 15×8 DataFrame
#  Row │ dna_pos  state   exon_number  exon_pos  intron_number  intron_pos  codon_number  codon_pos 
#      │ Int64    String  Int64?       Int64?    Int64?         Int64?      Int64?        Int64?    
# ─────┼────────────────────────────────────────────────────────────────────────────────────────────
#    1 │       1  UTR         missing   missing        missing     missing       missing    missing 
#    2 │       2  exon              1         1        missing     missing             1          1
#    3 │       3  exon              1         2        missing     missing             1          2
#    4 │       4  intron      missing   missing              1           1       missing    missing 
#    5 │       5  intron      missing   missing              1           2       missing    missing 
#    6 │       6  exon              2         3        missing     missing             1          3
#    7 │       7  exon              2         4        missing     missing             2          1
#    8 │       8  exon              2         5        missing     missing             2          2
#    9 │       9  intron      missing   missing              2           1       missing    missing 
#   10 │      10  exon              3         6        missing     missing             2          3
#   11 │      11  exon              3         7        missing     missing             3          1
#   12 │      12  exon              3         8        missing     missing             3          2
#   13 │      13  exon              3         9        missing     missing             3          3
#   14 │      14  UTR         missing   missing        missing     missing       missing    missing 
#   15 │      15  UTR         missing   missing        missing     missing       missing    missing 



@testset "dna_positions" begin
    @test dna_pos.dna_pos == collect(1:15)
    @test isequal(dna_pos.codon_number, [missing, 1, 1, missing, missing, 1, 2, 2, missing, 2, 3, 3, 3, missing, missing])
end
