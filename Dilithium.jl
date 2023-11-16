###  Dilithium
module Dilithium

using Nemo,Random

struct Param
    n::Int
    q::Int
    k::Int
    l::Int
    eta::Int
    R::zzModPolyRing       # R = ZZ/(q)ZZ
    mod::zzModPolyRingElem 
end

function Param(n::Int = 256, q::Int = 8380417, k::Int = 5, l::Int = 7, eta::Int = 2) 
    ispow2(n)                               || @error "n is not a power of 2"
    isprime(q)                              || @error "q is not prime"
    (1 == q % (2 * n))                      || @warn "NTT not supported by this choice of parameters"
    R, x = PolynomialRing(ResidueRing(ZZ, q), "x")
    mod = (gen(R) ^ n + 1)
    return Param(n, q, k, l, eta, R, mod)
end

struct Pk{T<:zzModPolyRingElem}
    A::Matrix{T}
    t::Matrix{T}
end

struct Sk{T<:zzModPolyRingElem}
    A::Matrix{T}
    t::Matrix{T}
    s1::Matrix{T}
    s2::Matrix{T}
end

function KeyGen(p::Param = Param())
    S,S1,S2 = MatrixSpace(p.R, p.k, p.l),MatrixSpace(p.R, p.l, 1),MatrixSpace(p.R, p.k, 1)

    A = Matrix(rand(S, 0:p.n-1))
    s1 = Matrix(rand(S1, -1:p.n, -p.eta:p.eta))
    s2 = Matrix(rand(S2, -1:p.n, -p.eta:p.eta))

    t = (A * s1 + s2) .% p.mod

    pk = Pk(A, t)
    sk = Sk(A, t, s1, s2)

    return (pk, sk)
end

"""
        array2ring(M::Matrix{zzModPolyRingElem}})

    returns  Matrix{Vector{Nemo.zzModRingElem}}, it converts ring element of M[i,j] into a vector c[i,j] of coefficients.
    Where c[i,j][1] belongs to the lowest degree coefficient and c[i,j][end] the highes degree coefficient.
"""
function ring2array(M::Matrix{zzModPolyRingElem})
    function elem2elem(e)
        #println(vcat(collect(coefficients(e)),zeros(base_ring(e),p.n-length(coefficients(e)))))
        return collect(coefficients(e))
    end 
    return elem2elem.(M)
    # returns type Matrix{Vector{Nemo.zzModRingElem}}
end 
"""
        pad_matrox(M::Matrix{Vector{Nemo.zzModRingElem}}, p::Param)

    returns again M::Matrix{Vector{Nemo.zzModRingElem}, but every coefficient vector is expanded to length(p.n)
"""
function pad_matrix(M::Matrix{Vector{Nemo.zzModRingElem}}, p::Param)
    function padvec(v)
        return vcat(v,zeros(parent(v[1]),p.n-length(v)))
    end 
    return padvec.(M)
    # returns type Matrix{Vector{Nemo.zzModRingElem}}
end 
"""
        array2ring(M::Matrix{Vector{Nemo.zzModRingElem}}, p::Param)

    returns  M::Matrix{zzModPolyRingElem}, it converts each coefficient vector to a ring element in p.R
"""
function array2ring(M::Matrix{Vector{Nemo.zzModRingElem}}, p::Param)
    function elem2elem(a)
        return sum(a.*[gen(p.R)^i for i = 0:length(a)-1])
    end 
    return elem2elem.(M)
    # returns type Matrix{zzModPolyRingElem}
end 
end # module Dilithium
