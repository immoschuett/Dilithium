include("Dilithium.jl")
include("Shake/shake.jl")
using .Dilithium, Test, Nemo, .SHAK3, Random, JSON
# gen challange: 

m = b"Hallo Welt" 
p = Dilithium.LV5
(pk,sk) = Dilithium.KeyGen(p)
sig = Dilithium.Sign(sk, m, p)
function export_challange(m ,p, pk, outfile = "c1.txt" )
    open(outfile, "w") do f
        println(f, "\n#PARAMETER: \nn = ", p.n, "\nq = ", p.q, "\nk = ", p.k, "\nl = ", p.l, "\neta = ", p.eta, " #η\ngamma1 = ", p.gamma1, " #γ_1\ngamma2 = ", p.gamma2, " #γ_2\ntau = ", p.tau, " #τ\nbeta = ", p.beta, " #β\nomega = ", p.omega, " #ω\nd = ", p.d)
        println(f, "\n#MESSAGE:")
        println(f, "\nem = ", json(m))
        println(f, "\n#PUBLIC_KEY:")
        println(f, "\nA = ", exportformat(pk.A,p))
        println(f, "\nt1 = ", exportformat(pk.t1,p))
        println(f, "\ntr = ", json(pk.tr))
        println(f, "\n#SIGNATURE: \n")
        println(f, "\nc0 = ", json(sig.c0))
        println(f, "\nz = ", exportformat(sig.z,p))
        println(f, "\nh = ", exportformat(sig.h,p))
    end # the file f is automatically closed after this block finishes
end
function exportformat(A,p=Dilithium.Param())
    if typeof(A)==Matrix{zzModPolyRingElem}
        A = Int.(lift.(Dilithium.ring2array(pk.A, p)))
    end
    S = size(A)
    B = Array{Array{Int}}(undef,S[1],S[2])
    for i=1:S[1],j=1:S[2]
        B[i,j] = A[i,j,:]
    end 
    B = json(B)
    if S[2]==1 
        return B[2:end-1]
    end 
    return B
end 
function importformat(In,p=Dilithium.Param(),nest=3)
    # input A as nested vector of Int
    if nest == 3 
        A = Array{Int}(undef,p.k,p.l,p.n)
        for i=1:p.k,j=1:p.l,k=1:p.n
            A[i,j,k] = In[i][j][k]
        end 
    end 
    if nest == 2
        A = Array{Int}(undef,length(In),p.n)
        for i=1:length(In),k=1:p.n
            A[i,k] = In[i][k]
        end 
    end 
    return A
end

@testset "In_export" begin
    m = b"Hallo Welt" 
    p = Dilithium.LV5
    (pk,sk) = Dilithium.KeyGen(p)
    sig = Dilithium.Sign(sk, m, p)
    @test  importformat(exportformat(pk.A))
end
export_challange(m ,p, pk)

# check that import(export) is the sameimportformat(exportformat(A))

m = b"Hallo Welt" 
p = Dilithium.LV5
(pk,sk) = Dilithium.KeyGen(p)
sig = Dilithium.Sign(sk, m, p)
importformat(exportformat(pk.A))