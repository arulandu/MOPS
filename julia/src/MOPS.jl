module MOPS

using SymPy

# Export functions
export parvalid, SFact, conjugate, arm, leg, Uhook, Lhook, GSFact, rho,
       subpar_check, partition_le, partition_lt, partition_ge, partition_gt, subPar, Par,
       JtoC, PtoC, JackIdentity, egen, evalJack, Jack, dominate, isdominate, Jack_c,
       monomial_sym, jack_sym, p_sym, c_sym, GBC, GBC_cont, GBC_cont_explicit, Jacobi_c, Jacobi

include("parvalid.jl")
include("sfact.jl")
include("conjugate.jl")
include("arm.jl")
include("leg.jl")
include("uhook.jl")
include("lhook.jl")
include("gsfact.jl")
include("rho.jl")
include("subpar_check.jl")
include("partition_compare.jl")
include("subpar.jl")
include("par.jl")
include("jtoC.jl")
include("ptoC.jl")
include("jack_identity.jl")
include("egen.jl")
include("evaljack.jl")
include("isdominate.jl")
include("dominate.jl")
include("jack_c.jl")
include("monomial_sym.jl")
include("gbc_cont.jl")
include("gbc_cont_explicit.jl")
include("gbc.jl")
include("jacobi_c.jl")
include("jacobi.jl")
include("jack.jl")

end # module MOPS
