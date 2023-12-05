include("Dilithium.jl")
include("Shake/shake.jl")
using .Dilithium, Test, Nemo, .SHAK3, Random, JSON
# gen challange: 

function manual_transpose(M,p)
    M = Int.(lift.(Dilithium.ring2array(M,p)))
    C = zeros(Int,p.l,p.k,p.n)
    for i = 1:p.l,j=1:p.k 
        C[i,j,:] = M[j,i,:]
    end 
    return C
end 
function export_challange(m ,p, pk, outfile = "c2.txt" )
    open(outfile, "w") do f
        println(f, "\n#PARAMETER: \nn = ", p.n, "\nq = ", p.q, "\nk = ", p.k, "\nl = ", p.l, "\neta = ", p.eta, " #η\ngamma1 = ", p.gamma1, " #γ_1\ngamma2 = ", p.gamma2, " #γ_2\ntau = ", p.tau, " #τ\nbeta = ", p.beta, " #β\nomega = ", p.omega, " #ω\nd = ", p.d)
        println(f, "\n#MESSAGE:")
        println(f, "\nem = ", json(m))
        println(f, "\n#PUBLIC_KEY:")
        #! this will transpose A since iterating over A will be columns !!! 
        println(f, "\nA = ", exportformat(manual_transpose(pk.A,p),p))
        println(f, "\nt1 = ", exportformat(pk.t1,p))
        println(f, "\ntr = ", json(pk.tr))
        println(f, "\n#SIGNATURE: \n")
        println(f, "\nc_tilde = ", json(sig.c0))
        println(f, "\nz = ", exportformat(sig.z,p))
        println(f, "\nh = ", exportformat(sig.h,p))
    end # the file f is automatically closed after this block finishes
end
function exportformat(In,p=Dilithium.Param())
    # ! do not change this !! 
    if typeof(In) == Matrix{Dilithium.T}
        In = Int.(lift.(Dilithium.ring2array(In,p)))
    end 
    S = size(In)
    B = Array{Array{Int}}(undef,S[1],S[2])
    for i=1:S[1],j=1:S[2]
        B[i,j] = In[i,j,:]
    end 
    B = json(B)
    if S[2]==1 
        return B[2:end-1]
    end 
    return B
end 
function importformat(In,p=Dilithium.Param(),nest=3)
    In = JSON.parse(In)
    # input A as nested vector of Int
    if nest == 3 
        A = Array{Int}(undef,p.k,p.l,p.n)
        for i=1:p.k,j=1:p.l,k=1:p.n
            #A[i,j,k] = In[j][i][k] # if own
            A[i,j,k] = In[i][j][k] # if Jens
        end 
    end 
    if nest == 2
        A = Array{Int}(undef,length(In),1,p.n)
        for i= 1:p.k,k=1:p.n
            A[i,1,k] = In[i][k]
        end 
    end 
    return A
end

m = b"Hallo Welt";
p = Dilithium.LV3
(pk,sk) = Dilithium.KeyGen(p);
pk.t1[1]
sig = Dilithium.Sign(sk, m, p);
sig.z
export_challange(m,p,pk)

@testset "In_export" begin
    m = b"Hallo Welt";
    p = Dilithium.LV3;
    (pk,sk) = Dilithium.KeyGen(p);
    sig = Dilithium.Sign(sk, m, p);
    @test importformat(exportformat(pk.A,p),p) == Dilithium.ring2array(pk.A,p)
    @test importformat(exportformat(pk.t1,p),p,2) == Dilithium.ring2array(pk.t1,p)
end
export_challange(m ,p, pk);

# check that import(export) is the sameimportformat(exportformat(A))


w = Dilithium.Vrfy(pk,m,sig,p)
www = Dilithium.bitpacking(w,p)
println(www)


m = b"Hallo Welt";
p = Dilithium.LV3;
(pk,sk) = Dilithium.KeyGen(p);
sig = Dilithium.Sign(sk, m, p);

exportformat(pk.A,p)
importformat(exportformat(pk.A,p),p) == Dilithium.ring2array(pk.A,p)
@test importformat(exportformat(pk.t1,p),p,2) == Dilithium.ring2array(pk.t1,p)