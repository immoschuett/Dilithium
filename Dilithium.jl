###  Dilithium
module Dilithium
include("Shake/shake.jl")
using Nemo, Random, .SHAK3

const AbstractBytes = Union{AbstractVector{UInt8},NTuple{N,UInt8} where N}

struct Param
    n::Int
    q::Int
    k::Int
    l::Int
    eta::Int
    gamma1::Int
    tau::Int
    beta::Int
    d::Int
    R::zzModPolyRing       # R = ZZ/(q)ZZ
    mod::zzModPolyRingElem
    seed::AbstractBytes
end

struct Pk{T<:zzModPolyRingElem}
    A::Matrix{T}
    t::Matrix{T}
end

struct Sk{T<:zzModPolyRingElem}
    # use 3-dimensional arrays 
    A::Matrix{T}
    t::Matrix{T}
    s1::Matrix{T}
    s2::Matrix{T}
end
# compressed PK, SK structs 
struct compPk{E<:zzModRingElem}
    A::Array{E,3}
    t1::Array{E,3}
    tr#::AbstractBytes
end

struct compSk{E<:zzModRingElem}
    A::Array{E,3}
    tr#::AbstractBytes
    s1::Array{E,3}
    s2::Array{E,3}
    t0::Array{E,3}
end

function Param(n::Int=256, q::Int=8380417, k::Int=5, l::Int=7, eta::Int=2, gamma1::Int = 2^17, tau::Int=39, beta::Int=78 , d::Int=13, seed::AbstractBytes=b"bads33d")
    ispow2(n)           || @error "n is not a power of 2"
    isprime(q)          || @error "q is not prime"
    (1 == q % (2 * n))  || @warn "NTT not supported by this choice of parameters"
    R, _ = PolynomialRing(ResidueRing(ZZ, q), "x")
    mod = (gen(R)^n + 1)
    return Param(n, q, k, l, eta, gamma1, tau, beta, d, R, mod, seed)
end

function KeyGen(p::Param=Param())
    BR = base_ring(p.R)

    A = rand(BR, p.k, p.l, p.n)
    s1 = BR.(rand(-p.eta:p.eta, p.l, 1, p.n))
    s2 = BR.(rand(-p.eta:p.eta, p.k, 1, p.n))

    A = array2ring(A, p)
    s1 = array2ring(s1, p)
    s2 = array2ring(s2, p)

    t = (A * s1 + s2) .% p.mod

    pk = Pk(A, t)
    sk = Sk(A, t, s1, s2)

    #return compress_key(sk,p)
    return (pk, sk)
end

function compress_key(sk::Sk, p::Param)
    # return compPk, compSk
    t1 = power2left.(ring2array(sk.t, p), Ref(p))
    t0 = power2right.(ring2array(sk.t, p), Ref(p))
    tr = SHAK3.shake256(vcat(p.seed, p.seed), UInt(32))
    return (compPk(ring2array(sk.A, p), t1, tr), compSk(ring2array(sk.A, p), tr, ring2array(sk.s1, p), ring2array(sk.s2, p), t0))
end



function Sign(csk::compSk, p::Param, m::AbstractBytes)

    BR = base_ring(p.R)
    A = array2ring(csk.A,p)
    s1 = array2ring(csk.s1,p)
    s2 = array2ring(csk.s2,p)
    mu = SHAK3.shake256(vcat(csk.tr, m), UInt(64))
    z, h = [], []
    while z == [] && h == []
        y = array2ring(BR.(rand(0:(2*p.gamma1 -1), p.l, 1, p.n)),p)
        w = A * y .% p.mod
        w = ring2array(w,p)
        w1 = highbits.(w, Ref(gamma1))
        # TODO reshape to Vector UInt8
        w1 = reinterpret(UInt8,reshape(w1,p.k*p.n))
        c = SHAK3.shake256(vcat(mu,w1),UInt(32))
        c = sampletoball(c,p)

        z = y+(s1.*c) .% p.mod 
        r0 = lowbits.((w - (s2.*c)) .% p.mod , Ref(gamma1))
        if (c_abs(z,p) >= (p.gamma1 - p.beta)) || (c_abs(r0,p) >= (p.gamma1 - p.beta)) 

        else 
            h = MakeHint
        end 
    end
end

function Vrfy()
end





"""
        array2ring(M::Matrix{zzModPolyRingElem}}, p::Param)

    returns  M::Array{Nemo.zzModRingElem},3}, it converts ring element of M[i,j] into a 3 dimensional Array C of coefficients.
    Where C[i,j,k] belongs to the k-th coefficient of the polynomial in M[i,j]. low k belong to low degree coefficients and high k to the highes degree coefficient.
"""
function ring2array(M::Matrix{zzModPolyRingElem}, p)
    (a, b) = size(M)
    padded = pad_matrix(collect.(coefficients.(M)), p)
    return [padded[i, j][k] for i = 1:a, j = 1:b, k = 1:p.n]
    # returns type Matrix{Vector{Nemo.zzModRingElem}}
end

function pad_matrix(M::Matrix{Vector{Nemo.zzModRingElem}}, p::Param)
    function padvec(v)
        return vcat(v, zeros(base_ring(p.R), p.n - length(v)))
    end
    return padvec.(M)
    # returns type Matrix{Vector{Nemo.zzModRingElem}}
end
"""
        array2ring(M::Array{Nemo.zzModRingElem},3}, p::Param)

    returns  M::Matrix{zzModPolyRingElem}, it converts the Array containing the coefficient vectors to a Matrix of ring elements in p.R
"""
function array2ring(M, p::Param)
    (a, b) = size(M)
    M2 = zeros(p.R, a, b)
    function elem2elem(a)
        return sum(a .* [gen(p.R)^i for i = 0:length(a)-1])
    end
    for i = 1:a, j = 1:b
        M2[i, j] = elem2elem(M[i, j, :])
    end
    return M2
    # returns type Matrix{zzModPolyRingElem}
end

function c_max(e::Z, p::Param) where {Z<:zzModRingElem}
    if lift(e) > div(p.q, 2)
        return abs(lift(e) - p.q)
    else
        return abs(lift(e))
    end
end
function c_abs(e,p)
    return c_max.(ring2array(e,p),Ref(p))
end 

# obsolet (splitted in left an right)
function power2rnd(r::Z, p::Param) where {Z<:zzModRingElem}
    mask = 2^p.d
    r0 = lift(r) % mask
    r = lift(r) % p.q
    # centered residue sys
    (r0 <= div(p.q, 2)) || (r0 -= mask)
    return (base_ring(p.R)((r - r0) >> p.d), base_ring(p.R)(r0))
end

function power2left(r::Z, p::Param) where {Z<:zzModRingElem}
    mask = 2^p.d
    r0 = lift(r) % mask
    r = lift(r) % p.q
    # centered residue sys
    (r0 <= div(p.q, 2)) || (r0 -= mask)
    return base_ring(p.R)((r - r0) >> p.d)
end

function power2right(r::Z, p::Param) where {Z<:zzModRingElem}
    mask = 2^p.d
    r0 = lift(r) % mask
    # centered residue sys
    (r0 <= div(p.q, 2)) || (r0 -= mask)
    return base_ring(p.R)(r0)
end

function decompose(r::Z, alpha, p::Param) where {Z<:zzModRingElem}
    r0 = lift(r) % alpha
    r = lift(r) % p.q
    # centered residue sys
    (r0 <= div(alpha, 2)) || (r0 -= alpha)
    if (r - r0) % p.q == p.q - 1
        r1 = 0
        r0 -= 1
    else
        r1 = divexact((r - r0) % p.q, alpha)
    end
    return (r1, r0)
end

#TODO move alpha into param !! 

function highbits(r::Z, alpha, p::Param) where {Z<:zzModRingElem}
    return decompose(r, alpha, p)[1]
end

function lowbits(r::Z, alpha, p::Param) where {Z<:zzModRingElem}
    return decompose(r, alpha, p)[2]
end

function makehint(z, r::Z, alpha, p) where {Z<:zzModRingElem}
    return highbits(r, alpha, p) != highbits(r + z, alpha, p)
end

function usehint(h, r::Z, alpha, p) where {Z<:zzModRingElem}
    p.q % alpha == 1 || @error "q != 1 mod alpha"
    m = divexact(p.q - 1, alpha)
    (r1, r0) = decompose(r, alpha, p)
    if h && r0 > 0
        return (r1 + 1) % m
    end
    if h && r0 <= 0
        return (r1 - 1 + m) % m
    end
    return r1
end

function sampletoball(rho::AbstractBytes, p::Param)
    hash = reinterpret(Int,rho) # ! change this in release
    c = vcat(zeros(Int,p.n-p.tau),rand(MersenneTwister(hash[1]),[-1,1],p.tau))
    shuffle!(MersenneTwister(hash[1]),c)
    # now we have a vector with  (p.n)  entries, tau +- 1
    c = base_ring(p.R)(c)
    # return polynomial of the vector 
    #(an element)
    return sum(c .* [gen(p.R)^i for i = 0:p.n-1])
end 



end # module Dilithium

