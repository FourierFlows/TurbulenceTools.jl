using TurbulenceTools.TwoDTurbTools

 n = 256
 ν = 1e-4
fi = 1.0
 μ = 1e0
nμ = -1
tf = 100

# 1/1e-2 = 100...
# 1/5e-3 = 200...
# 1/2e-3 = 500...

for ki = [16, 24]
  for dt = [2e-3, 1e-3]
    @printf("Running with dt=%3.0e, fi=%.2f, ki=%d, and μ=%4.0e...\n", 
      dt, fi, ki, μ)

    plotname = @sprintf("stochastic_fi%02d_n%d_ki%02d_mu%4.0e_dt%4.0e", 
      Int(fi*10), n, ki, μ, dt)

    runforcingproblem(fi=fi, ki=ki, tf=tf, dt=dt, withplot=true, 
      ns=5, n=n, ν=ν, μ=μ, nμ=nμ, stepper="FilteredRK4", plotname=plotname,
      withoutput=true, filename=plotname, stochastic=true)
  end
end
