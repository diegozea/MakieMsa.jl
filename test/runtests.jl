using Test


@testset "dna_location" begin
    @test parse_location("join(123..234,324..456)") == [[123,234],[324,456]]
    @test parse_location("1..200") == [[1,200]]    
end

dna_pos = dna_positions("2..3,6..8,10..13", [1,15])
@testset "dna_positions" begin
    @test dna_pos.dna_pos == collect(1:15)
    @test isequal(dna_pos.codon_number, [missing, 1, 1, missing, missing, 1, 2, 2, missing, 2, 3, 3, 3, missing, missing])
end
