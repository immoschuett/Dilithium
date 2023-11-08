###  Dilithium
using Nemo
include("Shake/shake.jl")
using Main.SHAK3

abstract type DILITHIUM_PARAMETER end
# Dilithium Parameter and Structs
struct D_Param<:DILITHIUM_PARAMETER # containing parameters
    n::Int
	q::Int
    k::Int
    l::Int
    eta::Int
end
struct D_Struct<:DILITHIUM_PARAMETER # containing parameters and structures
    P::D_Param                      # Parameter
    R::zzModPolyRing                # R = ZZ/(q)ZZ
    mod::zzModPolyRingElem          # x^n + 1 \in R
end 
# Dilithium Context, containing the seeds to build all elements. 
struct D_CTX<:DILITHIUM_PARAMETER
    rho::Array{UInt8, 1}
    rho_prime::Array{UInt8, 1}
    K::Array{UInt8, 1}
end
function init_Dctx(seed::Base.CodeUnits{UInt8, String})
    H = shake256(seed,0x000000080)
    return D_CTX(H[1:32],H[32:97],H[97:end])
end  
function init_param(n::Int=256,q=8380417::Int,k=10::Int,l=10::Int,eta=2::Int)
    ispow2(n) || @error "n is not a power of 2"
    isprime(q) || @error "q is not prime"
    (1 == q%(2*n)) || @warn "NTT not supported by this choice of parameters"
    return D_Param(n,q,k,l,eta)
end 
function init_ring(param::D_Param=init_param())
    R,x=PolynomialRing(ResidueRing(ZZ,param.q),"x")
    mod = (gen(R)^param.n + 1)
    return D_Struct(param,R,mod)
end 
function key_gen(DIL::D_Struct)
    param = DIL.P
    S = MatrixSpace(DIL.R,param.k,param.l)
    # samle according the given distributions
    S1 = MatrixSpace(DIL.R,param.l,1)
    S2 = MatrixSpace(DIL.R,param.k,1)   
    A = rand(S,0:param.n)  # ! degree >= 0, this excludes the 0 here ! 
    s1 = rand(S1,-1:param.n,-param.eta:param.eta)
    s2 = rand(S2,-1:param.n,-param.eta:param.eta)
    t = (A*s1+s2).% DIL.mod
    pk = (A,t)
    sk = (A,t,s1,s2)
    return (pk,sk)
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

#########################################################################################################################################
### EXAMPLES
orig_DiL = init_ring()
DIL = init_ring(init_param(7,6,10,10,2))
key = key_gen(DIL)
A = key[1][1]
digitalA = format(A,DIL)
#########################################################################################################################################
# shake 
shake128(b"ab",UInt(1000))
shake256(b"ab",UInt(912))
sha3_256(b"ab")
testseed = b"TeuleXOOIwXiRtofPsOFNbh2dHaVlAGZ"
init_Dctx(testseed)
#########################################################################################################################################
# expand from seed
expandA()

#########################################################################################################################################
#ntt
R = ResidueRing(ZZ,17)
P,x = PolynomialRing(R,"x")

zeta = R(2) # 8-th root of unity mod 17
testp1 = 3+2x+x^3
testp2 = 4+12x+x^2+7x^3
t1 = collect(coefficients(testp1))
T1 = vcat(t1,R.([0,0,0,0]))
t2 = collect(coefficients(testp2))
T2 = vcat(t2,R.([0,0,0,0]))
NTT(8,T1,zeta)

d = INTT(8,NTT(8,T1,zeta).*NTT(8,T2,zeta),zeta)
r =INTT(8,NTT(8,T1,zeta).*NTT(8,T2,zeta),zeta,true)
testp1*testp2
testp1*testp2.%(x^4+1)