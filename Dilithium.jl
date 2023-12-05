###  Dilithium
module Dilithium
include("Shake/shake.jl")
using Nemo, Random, .SHAK3

const AbstractBytes = Union{AbstractVector{UInt8},NTuple{N,UInt8} where N}
const Z = zzModRingElem
const T = zzModPolyRingElem

struct Param
    n::Int
    q::Int
    k::Int
    l::Int
    eta::Int                # η
    gamma1::Int             # γ_1
    gamma2::Int             # γ_2
    tau::Int                # τ
    beta::Int               # β
    omega::Int              # ω
    d::Int
    R::zzModPolyRing        # R = ZZ/(q)ZZ (x)
    mod::zzModPolyRingElem
    seed::AbstractBytes     # ρ'
end

struct PublicKey
    A::Matrix{T}
    t1::Matrix{T}
    tr::AbstractBytes
end

struct SecretKey
    A::Matrix{T}
    tr::AbstractBytes
    s1::Matrix{T}
    s2::Matrix{T}
    t0::Matrix{T}
end

struct Signature
    c0::Array{UInt8,1}
    z::Matrix{T}
    h::BitArray{3}
end

function Param(; n::Int=256, q::Int=8380417, k::Int=4, l::Int=4, eta::Int=2, gamma1::Int=2^17, gamma2::Int=divexact(8380416, 88), tau::Int=39, omega::Int=80, d::Int=13, seed::AbstractBytes=b"bads33d")::Param
    ispow2(n) || @error "n is not a power of 2"
    isprime(q) || @error "q is not prime"
    (1 == q % (2 * n)) || @warn "NTT not supported by this choice of parameters"
    tau <= n || @error "tau > n"
    R, _ = PolynomialRing(ResidueRing(ZZ, q), "x")
    mod = (gen(R)^n + 1)
    beta = 2 * tau
    return Param(n, q, k, l, eta, gamma1, gamma2, tau, beta, omega, d, R, mod, seed)
end

# NIST SECURITY LVL
const LV2 = Param()

const LV3 = Param(n=256, q=8380417, k=6, l=5, eta=4, gamma1=2^19, gamma2=divexact(8380416, 32), tau=49, omega=55)

const LV5 = Param(n=256, q=8380417, k=6, l=5, eta=2, gamma1=2^19, gamma2=divexact(8380416, 32), tau=60, omega=75)

function KeyGen(p::Param=Param())::Tuple{PublicKey,SecretKey}
    BR = base_ring(p.R)

    A = rand(BR, p.k, p.l, p.n)
    s1 = BR.(rand(-p.eta:p.eta, p.l, 1, p.n))
    s2 = BR.(rand(-p.eta:p.eta, p.k, 1, p.n))

    A = array2ring(A, p)
    s1 = array2ring(s1, p)
    s2 = array2ring(s2, p)

    t = (A * s1 + s2) .% p.mod

    t1 = power2left.(ring2array(t, p), Ref(p))
    t0 = array2ring(power2right.(ring2array(t, p), Ref(p)), p)
    tr = shake256(vcat(p.seed, reinterpret(UInt8, UInt32.(reshape(lift.(t1), p.k * p.n)))), UInt(32))

    t1 = array2ring(t1, p)
    return (PublicKey(A, t1, tr), SecretKey(A, tr, s1, s2, t0))
end

function Sign(sk::SecretKey, m::AbstractBytes, p::Param)::Signature
    BR = base_ring(p.R)
    mu = shake256(vcat(sk.tr, m), UInt(64))
    c0 = undef
    z, h = undef, undef
    while z == undef && h == undef
        # return y in S_gamma1 ^ l 
        y = array2ring(BR.(rand(-p.gamma1:p.gamma1, p.l, 1, p.n)), p)
        w = (sk.A * y) .% p.mod
        w1 = HighBits(w, p)
        c0 = shake256(vcat(mu, bitpacking(w1, p)), UInt(32))
        c = sampleinball(c0, p)
        z = (y + (sk.s1 .* c)) .% p.mod
        r0 = LowBits((w .- (sk.s2 .* c)) .% p.mod, p)
        # note Low,Higbits return ZZRingElem in centered representation
        if ((c_abs(z, p) >= (p.gamma1 - p.beta))) | ((c_abs_z(r0) >= (p.gamma2 - p.beta)))
            z, h = undef, undef
        else
            h = MakeHint((sk.t0 .* (c * (-1))) .% p.mod, (w - sk.s2 .* c + sk.t0 .* c) .% p.mod, p)
            if (c_abs((sk.t0 .* c) .% p.mod, p) >= p.gamma2) | (sum(h) > p.omega)
                z, h = undef, undef
            end
        end
    end
    return Signature(c0, z, h)
end

function Vrfy(pk::PublicKey, m::AbstractBytes, sig::Signature, p::Param)#::Bool
    mu = shake256(vcat(pk.tr, m), UInt(64))
    c = sampleinball(sig.c0, p)
    w1 = UseHint(sig.h, (pk.A * sig.z - ((pk.t1 .* c) .* (2^p.d))) .% p.mod, p)
    return (c_abs(sig.z, p) < (p.gamma1 - p.beta) &&
            (sum(sig.h) <= p.omega) &&
            (sig.c0 == shake256(vcat(mu, bitpacking(w1, p)), UInt(32))))
end

function bitpacking(w,p)::Array{UInt8,1}
    @assert ceil((p.q-1)/(2*p.gamma2)) == ((p.q-1)/(2*p.gamma2))
    pr = Int(ceil(log2(Int((p.q-1)/(2*p.gamma2))))) # packing radius
    # either 4 or 6 
    @assert (pr in [4,6])
    ret = UInt8[]
    if pr == 4
        for i = 1:p.k, k = 1:2:p.n
            a = UInt8(w[i,1,k])
            b = UInt8(w[i,1,k+1])
            c = (a << 4) | b
            push!(ret,c)
        end 
    end 
    if pr == 6 
        for i = 1:p.k, k = 1:3:p.n
            a = UInt8(w[i,1,k])
            b = UInt8(w[i,1,k+1])
            c = UInt8(w[i,1,k+2])
            d = (a << 2) | ( b >> 4)
            e = (b << 4) | c
            push!(ret,d,e)
        end 
    end 
    return ret
end


"""
        array2ring(M::Matrix{zzModPolyRingElem}}, p::Param)

    returns  M::Array{Nemo.zzModRingElem},3}, it converts ring element of M[i,j] into a 3 dimensional Array C of coefficients.
    Where C[i,j,k] belongs to the k-th coefficient of the polynomial in M[i,j]. low k belong to low degree coefficients and high k to the highes degree coefficient.
"""
function ring2array(M::Matrix{T}, p::Param)::Array{Z,3}
    (a, b) = size(M)
    padded = pad_matrix(collect.(coefficients.(M)), p)
    return [padded[i, j][k] for i = 1:a, j = 1:b, k = 1:p.n]
end

function pad_matrix(M::Matrix{Vector{Z}}, p::Param)::Matrix{Vector{Z}}
    function padvec(v)
        return vcat(v, zeros(base_ring(p.R), p.n - length(v)))
    end
    return padvec.(M)
end

function array2ring(M::Array{Z,3}, p::Param)::Matrix{T}
    (a, b) = size(M)
    M2 = zeros(p.R, a, b)
    function elem2elem(a)
        return sum(a .* [gen(p.R)^i for i = 0:length(a)-1])
    end
    for i = 1:a, j = 1:b
        M2[i, j] = elem2elem(M[i, j, :])
    end
    return M2
end
function array2ring(M::Array{Int,3}, p::Param)::Matrix{T}
    (a, b) = size(M)
    M2 = zeros(p.R, a, b)
    function elem2elem(a)
        return sum(a .* [gen(p.R)^i for i = 0:length(a)-1])
    end
    for i = 1:a, j = 1:b
        M2[i, j] = elem2elem(M[i, j, :])
    end
    return M2
end

# Functions for ZZRingElem
function c_abs_z(e::Array{ZZRingElem})
    return maximum(abs.(e))
end

# Functions for Z::zzModRingElem
function c_max(e::Z, p::Param)
    if lift(e) > div(p.q, 2)
        return abs(lift(e) - p.q)
    else
        return abs(lift(e))
    end
end
# obsolet (splitted in left an right)
function power2rnd(r::Z, p::Param)
    mask = 2^p.d
    r0 = lift(r) % mask
    r = lift(r) % p.q
    # centered residue sys
    (r0 <= div(p.q, 2)) || (r0 -= mask)
    return (base_ring(p.R)((r - r0) >> p.d), base_ring(p.R)(r0))
end
function power2left(r::Z, p::Param)
    mask = 2^p.d
    r0 = lift(r) % mask
    r = lift(r) % p.q
    # centered residue sys
    (r0 <= div(p.q, 2)) || (r0 -= mask)
    return base_ring(p.R)((r - r0) >> p.d)
end

function power2right(r::Z, p::Param)
    mask = 2^p.d
    r0 = lift(r) % mask
    # centered residue sys
    (r0 <= div(p.q, 2)) || (r0 -= mask)
    return base_ring(p.R)(r0)
end

function decompose(r::Z, alpha::Int, p::Param)
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

function highbits(r::Z, alpha::Int, p::Param)
    return decompose(r, alpha, p)[1]
end

function lowbits(r::Z, alpha::Int, p::Param)
    return decompose(r, alpha, p)[2]
end

function makehint(z::Z, r::Z, alpha::Int, p::Param)
    return highbits(r, alpha, p) != highbits(r + z, alpha, p)
end

function usehint(h, r::Z, alpha::Int, p::Param)
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

# Functions for T::zzModPolyRingElem
function MakeHint(z::Matrix{T}, r::Matrix{T}, p)
    z = ring2array(z, p)
    r = ring2array(r, p)
    return makehint.(z, r, Ref(2 * p.gamma2), Ref(p))
end

function UseHint(h, r::Matrix{T}, p)
    r = ring2array(r, p)
    return usehint.(h, r, Ref(2 * p.gamma2), Ref(p))
end

function HighBits(r::Matrix{T}, p)
    r = ring2array(r, p)
    return highbits.(r, Ref(2 * p.gamma2), Ref(p))
end

function LowBits(r::Matrix{T}, p)
    r = ring2array(r, p)
    return lowbits.(r, Ref(2 * p.gamma2), Ref(p))
end

function c_abs(e::Matrix{T}, p)
    return maximum(c_max.(ring2array(e, p), Ref(p)))
end

# additional functions
function sampleinball(rho::AbstractBytes, p::Param)
    ctx = SHAKEByteExtractor(SHAKE_256_CTX(), rho)

    rand_bytes = [extract_byte(ctx) for _ = 1:8]
    @assert(p.tau <= p.n)
    @assert(p.tau <= 64)
    sign_bits = vcat(digits.(rand_bytes, base=2, pad=8)...)[1:p.tau]
    c = zeros(Int, p.n)

    for i = p.n-p.tau:p.n-1
        j = extract_byte(ctx)
        while j > i
            j = extract_byte(ctx)
        end

        s = sign_bits[i-(p.n-p.tau)+1]

        c[i+1] = c[j+1]
        c[j+1] = (-1)^s
    end

    # ensure we generate only valid c with 
    @assert(sum(c .!= 0) == p.tau)

    return sum(c .* [gen(p.R)^i for i = 0:p.n-1])
end

end # module Dilithium
