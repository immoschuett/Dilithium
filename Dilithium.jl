###  Dilithium
using Nemo
include("Shake/shake.jl")
using Main.SHAK3

abstract type DILITHIUM_PARAMETER end
# Dilithium Parameter and Structs
struct D_Param<:DILITHIUM_PARAMETER
    n::Int
	q::Int
    k::Int
    l::Int
    eta::Int
end
struct D_Struct<:DILITHIUM_PARAMETER
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
    H = shake256(seed,128)
    return D_CTX(H[1:32],H[32:97],H[97:end])
end  
function init_param(n::Int,q::Int,k::Int,l::Int,eta::Int)
    ispow2(n) || @error "n is not a power of 2"
    isprime(q) || @error "q is not prime"
    (1 == q%(2*n)) || @warn "NNT not supported by this choice of parameters"
    return D_Param(n,q,k,l,eta)
end 
function init_ring(param::D_Param)
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
### EXAMPLES
DIL = init_ring(init_param(7,6,10,10,2))
A = key_gen(DIL)[1][1]
digitalA = format(A,DIL)



### OLD Workspace
#= Folded Distribution ?? 
function sample(fold,eta)
    # return a sample from the k-fold uniformly distr.
    sum(rand(-eta:eta,fold))
end
function sample_v(fold,eta,len)
    # return a sample array from the k-fold uniformly distr.
    sum(rand(-eta:eta,fold,len),dims=1)
end  
function convert()
# convert -ring elements to coefficient vectors w.r.t. our ordering.
#         -arrays of ring elements to corresponding arrays w.r.t. our ordering.
end =#

# Tests:
@time shake256(b"test",UInt(1000))
shake128(b"ab",UInt(1000))
sha3_256(b"ab")
sha3_512(b"ab124123")
testseed = b"TeuleXOOIwXiRtofPsOFNbh2dHaVlAGZ"
init_Dctx(testseed)
shake256(shake128(b"ab",UInt(1000)),UInt(10000))

length("TeuleXOOIwXiRtofPsOFNbh2dHaVlAGZTeuleXOOIwXiRtofPsOFNbh2dHaVlAGZTeuleXOOIwXiRtofPsOFNbh2dHaVlAGZTeuleXOOIwXiRtofPsOFNbh2dHaVlAGZTeuleXOOIwXiRtofPsOFNbh2dHaVlAGZTeuleXOOIwXiRtofPsOFNbh2dHaVlAGZTeuleXOOIwXiRtofPsOFNbh2dHaVlAGZTeuleXOOIwXiRtofPsOFNbh2dHaVlAGZTeuleXOOIwXiRtofPsOFNbh2dHaVlAGZ")
#TODO prealloc memory for shake! 


function expandA()

end 

function NTT(n::Int,f::Array{zzModRingElem, 1},zeta::zzModRingElem,Zeta=[zeta^i for i = 0:n]::Array{zzModRingElem, 1})
    # f is coefficient vector of polynomial
    ispow2(n) || @error "n not a power of 2"
    if n == 1
        return f
    end 
    #check zeta is n-th root of unity 
    # TODO 
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
    # shift and/or reduce with x^n/2 +1 
    if reduce
        return I[1:divexact(n,2)].-I[divexact(n,2)+1:end]
    end 
    return I
end 


using Hecke

#example NTT
using Nemo

R = ResidueRing(ZZ,17)

f = R.([2,1,0,7,0,0,0,0])
P,x = PolynomialRing(R,"x")

# achtung jetzt andere reihenfolge

eta = R(2) # 8-th root of unity
testp1 = 3+2x+x^3
testp2 = 4+12x+x^2+7x^3
t1 = collect(coefficients(testp1))
T1 = vcat(t1,R.([0,0,0,0]))
t2 = collect(coefficients(testp2))
T2 = vcat(t2,R.([0,0,0,0]))
NTT(8,T1,eta)

[testp1(eta^(i)) for i = 0:7]
t2 = collect(coefficients(testp2))
T2 = vcat(R.([0,0,0,0]),t2[end:-1:1])

eta = R(2) # 8-th primitive root mot 17
d = INTT(8,NTT(8,T1,eta).*NTT(8,T2,eta),eta,true)

testp1*testp2.%(x^4+1)
9-14 +17
t4 = vcat(collect(coefficients(testp1*testp2)),R.([0]))
T4 = t4[end:-1:1]
NTT(8,T4,eta)

2^7 % 17

INTT(8,NTT(8,T1,eta).*NTT(8,T2,eta),eta)

is_primitive(2,17)

inv(eta)
function MUL(A::NTT_Rep,B::NTT_Rep)

end 