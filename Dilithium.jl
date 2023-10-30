###  Dilithium
using Nemo
abstract type DILITHIUM_PARAMETER end
struct D_Param<:DILITHIUM_PARAMETER
    n::Int
	q::Int
    k::Int
    l::Int
    eta::Int
end
struct D_Struct<:DILITHIUM_PARAMETER
    P::D_Param
    R::zzModPolyRing
    mod::zzModPolyRingElem
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



