using TurbulenceTools.TwoDTurbTools, JLD2

include("sampleparams.jl")
 
@printf("Running stochastic forced turbulence for %d steps...\n", nt)

plotname = @sprintf("equilq_n%d_Ro%.1e_Rmu%.1e_Rnu_%.1e_ki%d", 
  n, Ro, Rμ, Rν, ki*L/2π)

prob, diags, output = initandrunproblem(L=L, n=n, fi=fi, ki=ki, tf=tf, dt=dt, 
  withplot=true, ns=ns, ν=ν, nν=nν, μ=μ, stepper="FilteredRK4", 
  plotname=plotname, withoutput=true, filename=plotname)
  
@printf "\n\nDone! max(Ro) = %.4f\n" maximum(abs.(prob.vars.q/f))
