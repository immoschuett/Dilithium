###  Dilithium
module Dilithium

using Nemo

struct Param
    n::Int
    q::Int
    k::Int
    l::Int
    eta::Int
    R::zzModPolyRing       # R = ZZ/(q)ZZ
    mod::zzModPolyRingElem 
end

function Param(n::Int = 256, q::Int = 8380417, k::Int = 10, l::Int = 10, eta::Int = 2) 
    ispow2(n) || @error "n is not a power of 2"
    isprime(q) || @error "q is not prime"
    (1 == q % (2 * n)) || @warn "NTT not supported by this choice of parameters"

    R, x = PolynomialRing(ResidueRing(ZZ, q), "x")
    mod = (gen(R) ^ n + 1)
    return Param(n, q, k, l, eta, R, mod)
end

struct Pk
    A::Matrix{zzModPolyRingElem}
    t::Matrix{zzModPolyRingElem}
end

struct Sk
    A::Matrix{zzModPolyRingElem}
    t::Matrix{zzModPolyRingElem}
    s1::Matrix{zzModPolyRingElem}
    s2::Matrix{zzModPolyRingElem}
end

function KeyGen(p::Param = Param())
    S = MatrixSpace(p.R, p.k, p.l)
    A = Matrix(rand(S, 0:p.n))
   
    S1 = MatrixSpace(p.R, p.l, 1)
    S2 = MatrixSpace(p.R, p.k, 1)
    s1 = Matrix(rand(S1, -1:p.n, -p.eta:p.eta))
    s2 = Matrix(rand(S2, -1:p.n, -p.eta:p.eta))

    t = (A * s1 + s2) .% p.mod

    pk = Pk(A, t)
    sk = Sk(A, t, s1, s2)

    return (pk, sk)
end
end # module Dilithium

module Temp
include("Shake/shake.jl")
using .SHAK3
using Nemo

# Dilithium Context, containing the seeds to build all elements. 
struct D_CTX
    rho::Array{UInt8, 1}
    rho_prime::Array{UInt8, 1}
    K::Array{UInt8, 1}
end

# small helper functions
function array2ring(a,R)
    x = gen(R)
    return sum(a.*[x^i for i = 0:length(a)-1])
end 
function ring2array(elem)
    return collect(coefficients(elem))
end 
function format(A,DIL)
    A = ring2array.(A)
    # filling up leading zeros
    for i in 1:size(A)[1],j in 1:size(A)[2]
        while length(A[i,j])<=DIL.P.n
            push!(A[i,j],zero(base_ring(DIL.R)))
        end 
    end 
    return A
end 

function init_Dctx(seed::Base.CodeUnits{UInt8, String})
    H = shake256(seed,0x000000080)
    return D_CTX(H[1:32],H[32:97],H[97:end])
end  

# TODO encript(sign), decrypt(vrfy), tests
#TODO prealloc memory for shake! Pull request to SHA.jl


# TODO:
function NTT(n::Int,f::Array{zzModRingElem, 1},zeta::zzModRingElem,Zeta=[zeta^i for i = 0:n]::Array{zzModRingElem, 1})
    # f is coefficient vector of polynomial
    ispow2(n) || @error "n not a power of 2"
    if n == 1
        return f
    end 
    #check zeta is n-th root of unity TODO
    # write f = g(x^2) + x*h(x)
    A = zeros(parent(f[1]),n)
    g = f[1:2:end]# even coeffs
    h = f[2:2:end]# odd coeffs
    B = NTT(divexact(n,2),g,Zeta[3],Zeta[1:2:end])  #B = NTT(divexact(n,2),g,zeta^2)
    C = NTT(divexact(n,2),h,Zeta[3],Zeta[1:2:end])  #C = NTT(divexact(n,2),h,zeta^2)
    for i = 1:divexact(n,2)
        A[i] = B[i] + Zeta[i]*C[i]                  #A[i] = B[i] + zeta^(i-1)*C[i]
        A[divexact(n,2)+i] = B[i] - Zeta[i]*C[i]    #A[divexact(n,2)+i] = B[i] - zeta^(i-1)*C[i]
    end 
    return A
end 

function INTT(n,f,zeta,reduce=false)
    I = inv(parent(f[1])(n)).*NTT(n,f,inv(zeta))
    # shift and/or reduce with x^n/2 +1 i.e subtract the higher registers from the lower registers, since  x^n \cong -1 
    if reduce  
        return I[1:divexact(n,2)].-I[divexact(n,2)+1:end]
    end 
    return I
end 

function expandA(D=init_Dctx(testseed)::D_CTX,DIL=init_ring(init_param(256,8380417,10,10,2))::D_Struct)
    # ! WARN this is not identical with the reference (different endianess / mapping into Zq, PROTOTYPE ONLY
    S = MatrixSpace(DIL.R,DIL.P.k,DIL.P.l)
    A = format(zero(S),DIL)
    rho = D.rho
    push!(rho,0x00)
    push!(rho,0x00)
    # compute a_ij 
    n = 256 # hardcode for now, later from D_param
    q = 8380417
    byte_number = Int(ceil(log2(q))) # number of bytes needed per coeff in NTT_REP
    for i in 1:DIL.P.k, j in 1:DIL.P.l
        b = 256*i+j 
        # watch out for byte order. (later) little endian order
        rho[32] = UInt8(b & 0x00ff )
        rho[33] = UInt8((b & 0xff00) >> 8)
        extract = shake128(rho,UInt(byte_number*n)) # per polynomial
        # reinterprete three bytes as value in Z/qZ
        # i.e reinterprete as integer and reduce mod q (later we can change this)
        for k = 0:n-1
            A[i,j][k+1]  = base_ring(DIL.R)(extract[3*k+1]^(extract[3*k+2]<<8)^(extract[3*k+3]<<16))
        end 
    end 
    return A
end 

end # module test
