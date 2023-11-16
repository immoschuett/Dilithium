include("Dilithium.jl")
using .Dilithium, Test, Nemo

@testset "conversion" begin
    for i = 1:10
        p = Dilithium.Param(16, 7, 5, 7,2) 
        A = Dilithium.KeyGen(p)[1].A
        @test Dilithium.array2ring(Dilithium.ring2array(A),p) == A
    end 
end

