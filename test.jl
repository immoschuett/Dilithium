include("Dilithium.jl")
include("Shake/shake.jl")
using .Dilithium, Test, Nemo, .SHAK3, Random


# TODO adjust hashing
# ! TODO LIST longterm
# TODO use NTT, 
# TODO use Pointwise mult. for 2+1 dim. Array with NTT 
# TODO compare with known answer testset
# TODO check secure randomness 

# Workspace 
#=
p = Dilithium.LV5
ck = Dilithium.compress_key( Dilithium.KeyGen(p)[2], p);
sg = Dilithium.Sign(ck[2], p, b"hallo");
typeof(sg[1])
typeof(sg[2])
typeof(sg[3])
println(sg[3])

sg[]
Dilithium.Vrfy(ck[1], b"hallo",sg, p)
=#
# Tesspace
@testset "verfysign" begin 
    for i = 1:3
        p = Dilithium.LV5
        (pk, sk) = Dilithium.KeyGen(p);
        sig = Dilithium.Sign(sk, b"hallo", p);
        @test Dilithium.Vrfy(pk, b"hallo",sig, p) == true
        sig.c0[1] += 1
        @test Dilithium.Vrfy(pk, b"hallo",sig, p) == false
    end 
    for i = 1:3
        p = Dilithium.LV3
        (pk, sk) = Dilithium.KeyGen(p);
        sig = Dilithium.Sign(sk, b"hallo", p);
        @test Dilithium.Vrfy(pk, b"hallo",sig, p) == true
        sig.c0[1] += 1
        @test Dilithium.Vrfy(pk, b"hallo",sig, p) == false
    end 
end 

@testset "hints" begin
    # check ID UseHintq(MakeHintq(z, r, α), r, α) = HighBitsq(z+r,alpha) 
    # if ||z||<= alpha/2
    #A = Dilithium.KeyGen(p)[1].A+ z, alpha).
    for i = 1:1000
        # q = 8380417
        # q > 2α
        # 1 = q%α
        alphahalf = gamma2 = 261888
        alpha = gamma1 = 2 * gamma2
        #alpha = 4092
        #alphahalf = 2046
        p = Dilithium.Param()
        z = base_ring(p.R)(rand(Int) % alphahalf) # then abs(z) <= alphahalf
        r = base_ring(p.R)(rand(Int) % alpha)
        @test Dilithium.usehint(Dilithium.makehint(z, r, alpha, p), r, alpha, p) == Dilithium.highbits(r + z, alpha, p)
        if Dilithium.usehint(Dilithium.makehint(z, r, alpha, p), r, alpha, p) != Dilithium.highbits(r + z, alpha, p)
            println("failes:")
            println(z)
            println(r)
        end
    end
end


@testset "SHAKE" begin
    # test some official testvectors from https://csrc.nist.gov/Projects/Cryptographic-Algorithm-Validation-Program/Secure-Hashing
    @test SHAK3.shake128(b"", UInt(16)) == hex2bytes("7f9c2ba4e88f827d616045507605853e")
    @test SHAK3.shake128(codeunits("0"^167), UInt(32)) == hex2bytes("ff60b0516fb8a3d4032900976e98b5595f57e9d4a88a0e37f7cc5adfa3c47da2")

    @test SHAK3.shake256(b"", UInt(32)) == hex2bytes("46b9dd2b0ba88d13233b3feb743eeb243fcd52ea62b81b82b50c27646ed5762f")
    @test SHAK3.shake256(codeunits("0"^135), UInt(32)) == hex2bytes("ab11f61b5085a108a58670a66738ea7a8d8ce23b7c57d64de83eaafb10923cf8")

end

@testset "conversion" begin
    for i = 1:10
        p = Dilithium.Param()
        A = Dilithium.KeyGen(p)[1].A
        @test Dilithium.array2ring(Dilithium.ring2array(A, p), p) == A
        @test Dilithium.array2ring(Dilithium.ring2array(zeros(p.R, 5, 7), p), p) == zeros(p.R, 5, 7)
    end
    # why error here?
end
