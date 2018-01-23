using TurbulenceTools.TwoDTurbTools.StochasticForcingProblems
using JLD2

 L = 1600e3
 n = 256

Ro = 0.1      # Rossby number
Rμ = 1e-4     # Ratio of drag and inertial time-scale
Rν = 1.0      # Ratio of viscous and inertial time-scale
 f = 1e-4     # Planetary vorticity

tq = 1/(f*Ro) # Eddy turnover time-scale
tμ = 1/(f*Rμ) # Bottom drag time-scale
tν = 1/(f*Rν) # Viscous time-scale

 μ = 1/tμ
kν = 2π/L * n/3 
 ν = 1/(kν*tν)
ki = 16*2π/L
fi = f*Ro/ki * sqrt(μ) # P = fi²

dt = tq/5
tf = 10tμ
nt = round(Int, tf/dt)
ns = 10

@printf("Running stochastic forced turbulence for %d steps...\n", nt)

plotname = @sprintf("equilq_n%d_Ro%.1e_Rmu%.1e_ki%d", n, Ro, Rμ, ki*L/2π)
  
prob, diags, filename = initandrunproblem(L=L, n=n, fi=fi, ki=ki, tf=tf, dt=dt, 
  withplot=true, ns=ns, ν=ν, μ=μ, stepper="FilteredRK4", plotname=plotname,
  withoutput=true, filename=plotname)

@printf "\n\nDone! max(Ro) = %.4f\n" maximum(abs.(prob.vars.q/f))

jldopen(filename, "r+") do file
  file["forcingparams/fi"] = fi
  file["forcingparams/ki"] = ki
end
