###  Dilithium
module Dilithium
export Param,Pk,Sk,KeyGen
include("Shake/shake.jl")
using .SHAK3,Nemo

struct D_CTX{T <: UInt8}
    rho::Array{T, 1}
    rho_prime::Array{T, 1}
    K::Array{T, 1}
    tr::Array{T, 1}
end

struct Param
    n::Int
    q::Int
    k::Int
    l::Int
    d::Int
    eta::Int
    R::zzModPolyRing       # R = ZZ/(q)ZZ
    mod::zzModPolyRingElem 
    ctx::D_CTX
end

struct Pk{T <: zzModPolyRingElem}  # later Matrix{Vector{zzModRingElem}}
    A:: Matrix{T}
    t::Vector{T}#debug only
    t1::Vector{T}
end

struct Sk{T <: zzModPolyRingElem}
    A:: Matrix{T}#also possible to use seed rho here
    t::Vector{T} #debug only
    t0::Vector{T}
    s1::Vector{T}
    s2::Vector{T}
end

#core functions

function init_ctx(seed::Base.CodeUnits{UInt8, String})
    H = shake256(seed,0x0000000a0)
    return D_CTX(H[1:32],H[32:97],H[97:128],H[129:end])
end  

function Param(n::Int = 256, q::Int = 8380417, k::Int = 5, l::Int = 7, d::Int = 13, eta::Int = 2, seed=b"") 
    ispow2(n)               || @error "n is not a power of 2"
    isprime(q)              || @error "q is not prime"
    (1 == q%(2*n))          || @warn  "NTT not supported by this choice of parameters"
    k < 256                 || @error "dimension k too high for sampling algorithm"
    l < 256                 || @error "dimension l too high for sampling algorithm"
    eta in [2,4]            || @error "eta has to be 2 or 4"
    R, _ = PolynomialRing(ResidueRing(ZZ, q), "x")
    mod = (gen(R) ^ n + 1)
    return Param(n, q, k, l, d, eta, R, mod, init_ctx(seed))
end

function KeyGen(p::Param = Param())
    A = expandA(p)
    s1,s2 = expandS(p::Param)
    t = (A * s1 + s2) .% p.mod
    pk = Pk(A,t, power2rnd_array(t,p)[1])
    sk = Sk(A, t, power2rnd_array(t,p)[2], s1, s2)
    return (pk, sk)
end

function Sign(sk,m,p)
    u = shake256(vcat(p.ctx.tr,m),0x0000000000000040)
    return u
end 

# small helper functions
# TODO make the helper functions nice and efficient
function array2ring(a,R)
    # return Matrix{zzModPolyRingElem} in the given ring R
    function a2r(a)
        return sum(a.*[gen(R)^i for i = 0:length(a)-1])
    end 
    return a2r.(a)
end 
function ring2array(elem)
    #return vector if coefficients
    return collect(coefficients(elem))
end 
function format(A,p::Param)# where T<:Union{Matrix{Vector{zzModRingElem}},Matrix{zzModPolyRingElem}}
    A = ring2array.(A)
    # filling up leading zeros
    for i in 1:size(A)[1],j in 1:size(A)[2]
        while length(A[i,j])<=p.n # ! ineff
            push!(A[i,j],zero(base_ring(p.R)))
        end 
    end 
    return A
end 

function power2rnd_array(a,p)
    b = deepcopy(a)
    c = deepcopy(a)
    for i=1:length(a)
        q,r = divrem(p.R(a[i]),p.R(2^p.d))
        b[i] = q
        c[i] = r
    end 
    return b,c
end 


function expandA(p::Param)
    S = MatrixSpace(p.R, p.k, p.l)
    A = format(zero(S),p) #k*l matrix with 
    rho = p.ctx.rho
    push!(rho,0x00,0x00)
    # compute a_ij 
    byte_number = Int(ceil(log2(p.q)/8)) # number of bytes needed per coeff in NTT_REP
    # evaluate a_ij
    for i in 1:p.k, j in 1:p.l  
        # append 2 nonce bytes = 256*i+j 
        # i,j < 256 by parameter choice
        rho[32] = UInt8(i)
        rho[33] = UInt8(j)
        # absorb rho | b into shake
        extract = shake128(rho,UInt(byte_number*p.n))
        # reinterprete byte_number bytes as value in ZZ/(q)ZZ
        for k = 0:p.n-1
            s = parse(UInt,"0x"*bytes2hex(extract[k*byte_number+1:k*byte_number+byte_number]))
            s &= 2^(8*byte_number) -1
            A[i,j][k+1] = base_ring(p.R)(s)
        end 
    end 
    return array2ring(A,p.R)
end 

function expandS(p::Param)
    S = MatrixSpace(p.R, p.k, p.l)
    s1 = format(zero(S),p)[1,:] # in R^l
    s2 = format(zero(S),p)[:,1] # in R^k

    rho_prime = p.ctx.rho_prime
    push!(rho_prime,0x00,0x00)
    p.l+p.k < 2^8 || @error "not implemented l+k does not fit into one byte"
    if p.eta ==2 # 1/2 byte per coefficient is needed
        # TODO
        for i=1:p.l
            # fill s1 
            rho_prime[end] = UInt8(i) # ? little endian 
            extract = shake256(rho_prime,UInt(divexact(p.n,2)))
            for j in 1:2:p.n
                s1[i][j]    = base_ring(p.R)(p.eta - (extract[divexact(j+1,2)] & 0xf)%5)
                s1[i][j+1]  = base_ring(p.R)(p.eta - ((extract[divexact(j+1,2)] >> 4) & 0xf)%5)
            end 
        end 
        # fill s2 
        for i=1:p.k
            rho_prime[66] = UInt8(p.l+i)
            extract = shake256(rho_prime,UInt(divexact(p.n,2)))
            for j in 1:2:p.n
                s2[i][j] = base_ring(p.R)(p.eta - (extract[divexact(j+1,2)] & 0xf)%5)
                s2[i][j+1] = base_ring(p.R)(p.eta - ((extract[divexact(j+1,2)] >> 4) & 0xf)%5)
            end 
        end 
    end 
    return array2ring(s1,p.R),array2ring(s2,p.R)
end 

function expandY()
    # TODO
end 

end # module Dilithium
