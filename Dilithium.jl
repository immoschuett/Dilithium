###  Dilithium
module Dilithium
include("Shake/shake.jl")
using Nemo,Random,.SHAK3

const AbstractBytes = Union{AbstractVector{UInt8},NTuple{N,UInt8} where N}

struct Param
    n::Int
    q::Int
    k::Int
    l::Int
    eta::Int
    R::zzModPolyRing       # R = ZZ/(q)ZZ
    mod::zzModPolyRingElem 
    seed::AbstractBytes
end

# seeds to produce A,s1,s2 from shake256 
struct Seed
    a#::AbstractBytes
    s#::AbstractBytes
    # tr :: TODO
end 

function Param(n::Int = 256, q::Int = 8380417, k::Int = 5, l::Int = 7, eta::Int = 2, seed::AbstractBytes = b"bads33d") 
    ispow2(n)                               || @error "n is not a power of 2"
    isprime(q)                              || @error "q is not prime"
    (1 == q % (2 * n))                      || @warn "NTT not supported by this choice of parameters"
    R, x = PolynomialRing(ResidueRing(ZZ, q), "x")
    mod = (gen(R) ^ n + 1)
    return Param(n, q, k, l, eta, R, mod, seed)
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

    S = subseeds(p.seed)
    A = expandA(p,S)
    s1,s2 = expandS(p,S)

    t = (A * s1 + s2) .% p.mod

    pk = Pk(A, t)
    sk = Sk(A, t, s1, s2)

    return (pk, sk)
end

function pow2rnd(t,d)

end 

"""
        array2ring(M::Matrix{zzModPolyRingElem}})

    returns  M::Matrix{Vector{Nemo.zzModRingElem}}, it converts ring element of M[i,j] into a vector c[i,j] of coefficients.
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
        return vcat(v,zeros(base_ring(p.R),p.n-length(v)))
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

function subseeds(seed::AbstractBytes)
    H = SHAK3.shake256(seed,0x000000080)
    return Seed(H[1:64],H[65:end])
end 


function expandA(p::Param, s::Seed)
    A = pad_matrix(ring2array(zeros(p.R,p.k,p.l)),p)
    byte_number = Int(ceil(log2(p.q)/8)) # exact number of bytes needed
    #Types =         [UInt8, UInt16, UInt32, UInt64, UInt128]
    #byte=number        1       2      3-4      5-8    9-16
    #log2(byte_number)  0       1      (1 2]   (2,3]   (3 4]
    Types =         [UInt8, UInt16, UInt32,  UInt32, UInt64, UInt64, UInt64, UInt64, UInt128, UInt128, UInt128, UInt128, UInt128, UInt128, UInt128, UInt128]
    Typ = Types[byte_number]
    mask = Typ(2)^(8*byte_number)-1  # 2^bitnumber -1
    extract = reinterpret(Typ,shake256(s.a,UInt(16*byte_number*p.n*p.k*p.l))) # estimated number of bytes needed for rejection sampling
    read = 1
    ptr = pointer(extract)
    len = length(extract)
    # set the coefficient vector in A[i,j]
    for i in 1:p.k, j in 1:p.l , k in 1:p.n
        # rejection sample byte_pnumber bytes from extract 
        # reinterpret, mask with 2^bitnumber-1  
        while (unsafe_load(ptr,read)&mask) > p.q
            read+=1
            read < len || @error "buffer overread"
        end 
        A[i,j][k] = base_ring(p.R)(unsafe_load(ptr,read)&mask) 
        read+=1
    end 
    return array2ring(A,p)
end

function expandS(p::Param, s::Seed)
    s1 = pad_matrix(ring2array(zeros(p.R,p.l,1)),p)
    s2 = pad_matrix(ring2array(zeros(p.R,p.k,1)),p)

    byte_number = Int(ceil(log2(p.eta)/8)) # exact number of bytes needed
    #Types =         [UInt8, UInt16, UInt32, UInt64, UInt128]
    #byte=number        1       2      3-4      5-8    9-16
    #log2(byte_number)  0       1      (1 2]   (2,3]   (3 4]
    Types =         [UInt8, UInt16, UInt32,  UInt32, UInt64, UInt64, UInt64, UInt64, UInt128, UInt128, UInt128, UInt128, UInt128, UInt128, UInt128, UInt128]
    Typ = Types[byte_number]
    mask = Typ(2)^(Int(ceil(log2(2*p.eta+1))))-1  # bits needed to  2*η 
    extract = reinterpret(Typ,shake128(s.s,UInt(16*byte_number*p.n*(p.k+p.l)))) # estimated number of bytes needed for rejection sampling
    read = 1
    ptr = pointer(extract)
    len = length(extract)
    # fill s1
    for i in 1:p.l, m in 1:p.n
        # rejection sample byte_pnumber bytes from extract 
        # reinterpret, mask with 2^bitnumber-1  
        while (unsafe_load(ptr,read)&mask) > 2* p.eta  # so we know we sampple values < 2 η  + 1
            read+=1
            read < len || @error "buffer overread"
        end 
        s1[i][m] = base_ring(p.R)(p.eta - unsafe_load(ptr,read)&mask) 
        read+=1
    end 
    for i in 1:p.k, m in 1:p.n
        # rejection sample byte_pnumber bytes from extract 
        # reinterpret, mask with 2^bitnumber-1  
        while (unsafe_load(ptr,read)&mask) >  2* p.eta
            read+=1
            read < len || @error "buffer overread"
        end 
        s2[i][m] = base_ring(p.R)(p.eta - unsafe_load(ptr,read)&mask) 
        read+=1
    end 
    return (array2ring(s1,p),array2ring(s2,p))
end


end # module Dilithium
