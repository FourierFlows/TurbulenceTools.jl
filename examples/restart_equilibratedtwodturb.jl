using TurbulenceTools.TwoDTurbTools, JLD2

include("sampleparams.jl")
 
filename = joinpath("data", 
  @sprintf("equilq_n%d_Ro%.1e_Rmu%.1e_Rnu_%.1e_ki%d", n, Ro, Rμ, Rν, ki*L/2π))
  
prob, diags, output = restartchanproblem(filename; tf=tf, ns=ns,
  plotname=filename, withplot=true)
  
@printf "\n\nDone! max(Ro) = %.4f\n" maximum(abs.(prob.vars.q/f))
