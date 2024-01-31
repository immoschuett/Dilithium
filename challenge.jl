include("Dilithium.jl")
include("Shake/shake.jl")
using .Dilithium, Test, Nemo, .SHAK3, Random, JSON
# gen challenge: 

function manual_transpose(M,p)
    M = Int.(lift.(Dilithium.ring2array(M,p)))
    C = zeros(Int,p.l,p.k,p.n)
    for i = 1:p.l,j=1:p.k 
        C[i,j,:] = M[j,i,:]
    end 
    return C
end 

function export_challenge(m ,p, pk, sig, outfile = "c1.txt" )
    open(outfile, "w") do f
        println(f, "\n#PARAMETER: \nn = ", p.n, "\nq = ", p.q, "\nk = ", p.k, "\nl = ", p.l, "\neta = ", p.eta, " #η\ngamma1 = ", p.gamma1, " #γ_1\ngamma2 = ", p.gamma2, " #γ_2\ntau = ", p.tau, " #τ\nbeta = ", p.beta, " #β\nomega = ", p.omega, " #ω\nd = ", p.d)
        println(f, "\n#MESSAGE:")
        println(f, "\nem = ", json(m))
        println(f, "\n#PUBLIC_KEY:")
        # this will transpose A since iterating over A will be columns !!! 
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
    # do not change this: 
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
    # WARNING  do not use this function!
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
        d = length(In)
        A = Array{Int}(undef,d,1,p.n)
        for i= 1:d,k=1:p.n
            A[i,1,k] = In[i][k]
        end 
    end 
    return A
end

#=
# PARSER:

# parser: / copy this to a c3.jl file.
p = Dilithium.LV3
A = importformat(string(A),p,3)
em = UInt8.(em)
c_tilde = UInt8.(c_tilde)
t1 = importformat(string(t1),p,2)
tr = UInt8.(tr)
z =  importformat(string(z),p,2)
h = Bool.(importformat(string(h),p,2))
z
# type pk: 
pk = Dilithium.PublicKey(Dilithium.array2ring(A,p),Dilithium.array2ring(t1,p),tr)
# type signature:
sig = Dilithium.Signature(c_tilde, Dilithium.array2ring(z,p), h)

# verify:
Dilithium.Vrfy(pk, em, sig, p)

#

# gen challenge:
m = b"Hallo Welt";
p = Dilithium.LV3
(pk,sk) = Dilithium.KeyGen(p);
pk.t1[1]
sig = Dilithium.Sign(sk, m, p);
sig.z
export_challenge(m,p,pk)

=#

# example generation:
# m = b"'Hearing voices no one else can hear isn’t a good sign, even in the wizarding world.' — Ron Weasley"
# #! IMPORTANT: later change to utf-8
# p = Dilithium.Param(n=16, k=23, l=23, tau=12)
# (pk, sk) = Dilithium.KeyGen(p);
# pk.t1[1]
# sig = Dilithium.Sign(sk, m, p);
# Dilithium.Vrfy(pk, m, sig, p) == true
# sig.z
# export_challenge(m, p, pk, "n=16.txt")

# All taus from moodle
taus = [4,5,6,7,8,9,10]
words_b32 = shuffle([])

for (idx,tau) in enumerate(taus)
    # Choose random 5-7 letter word, change in some way
    m = words_b32[idx]

    param = Dilithium.Param(tau = tau)
    (pk, sk) = Dilithium.KeyGen(param)
    sig = Dilithium.Sign(sk, m, param)

    Dilithium.Vrfy(pk, m, sig, param) == true || @error "Could not vrfy signature!"

    file_name = replace("challenge_tau_{}.txt", "{}" => string(tau, base = 10, pad = 2))
    export_challenge(m, param, pk, sig, file_name)
end

struct TestParam
    n :: Int
    k :: Int
    l :: Int
end

# More messages
messages = shuffle([])

# All params from moodle
parameters = [
    TestParam(16,18,18),
    TestParam(32,9,6),
    TestParam(16,19,19),
    TestParam(32,9,5),
    TestParam(64,5,5),
    TestParam(16,20,20),
    TestParam(16,21,21),
    TestParam(16,22,22),
    TestParam(16,23,23),
    TestParam(64,6,6),
    TestParam(64,6,5),
    TestParam(16,24,24),
    # TestParam(n, k, l),
]

m_idx = 1
for (; n, k, l) = parameters
    tau = Dilithium.LV2.tau
    if tau > n 
        tau = divexact(3 * n, 4)
    end

    m0 = messages[m_idx]

    # swap first and last character
    m = m0[end] * m[2:end-1] * m0[1]

    param = Dilithium.Param(n = n, k = k, l = l, tau = tau)
    (pk, sk) = Dilithium.KeyGen(param)
    sig = Dilithium.Sign(sk, m, param)

    Dilithium.Vrfy(pk, m, sig, param) == true || @error "Could not vrfy signature!"

    file_name = replace("challenge_uf_{}.txt", "{}" => string(m_idx, base = 10, pad = 2))
    export_challenge(m, param, pk, sig, file_name)

    global m_idx = m_idx + 1
end
