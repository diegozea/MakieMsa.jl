using Test

testset "dna_location" begin
    @test exons_from_location("join(123..234,324..456)") == [[123,234],[324,456]]
end
