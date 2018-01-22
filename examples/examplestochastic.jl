using TurbulenceTools.TwoDTurbTools

 n = 256
 ν = 1e-4
fi = 1.0
 μ = 1e0
tf = 100
ki = 24
dt = 2e-3

@printf("Running with dt=%3.0e, fi=%.2f, ki=%d, and μ=%4.0e...\n", 
  dt, fi, ki, μ)

plotname="example"

runforcingproblem(fi=fi, ki=ki, tf=tf, dt=dt, withplot=true, 
  ns=5, n=n, ν=ν, μ=μ, stepper="FilteredRK4", plotname=plotname,
  withoutput=true, filename=plotname, stochastic=true)
