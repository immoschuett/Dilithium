module test 
using Nemo
# Dilithium Context, containing the seeds to build all elements. 
struct D_CTX
    rho::Array{UInt8, 1}
    rho_prime::Array{UInt8, 1}
    K::Array{UInt8, 1}
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
end # module test
