include("Dilithium.jl")
include("Shake/shake.jl")
using .Dilithium, Test, Nemo, .SHAK3, Random
# gen challange: 

m = b"hallo"
p = Dilithium.LV5
(pk, sk) = Dilithium.KeyGen(p)
sig = Dilithium.Sign(sk, m, p)
#Dilithium.Vrfy(pk, bm, sig, p)



outfile = "c1.txt"
open(outfile, "w") do f
    println(f, "\nPARAMETER: \nn = ", p.n, "\nq = ", p.q, "\nk = ", p.k, "\nl = ", p.l, "\nη = ", p.eta, "\nγ_1 = ", p.gamma1, "\nγ_2 = ", p.gamma2, "\nτ = ", p.tau, "\nβ = ", p.beta, "\nω = ", p.omega, "\nd = ", p.d)
    println(f, "\nMESSAGE: \n")
    println(f, "m = ", m)
    println(f, "\nPUBLIC_KEY:")
    println(f, "\nA = \n", Int.(lift.((Dilithium.ring2array(pk.A, p)))))
    println(f, "\nt1 = \n", Int.(lift.(Dilithium.ring2array(pk.t1, p))))
    println(f, "\ntr = \n", pk.tr)
    println(f, "\nSIGNATURE: \n")
    println(f, "\nc0 = \n", sig.c0)
    println(f, "\nz = \n", Int.(lift.(Dilithium.ring2array(sig.z, p))))
    println(f, "\nh = \n", sig.h)
end # the file f is automatically closed after this block finishes