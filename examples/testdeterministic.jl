using TurbulenceTools.TwoDTurbTools

 n = 512
 ν = 1.0e-3
fi = 1
ki = 16
 μ = 0.01 #0.1*(2π)^2
tf = 10


# 1/1e-2 = 100...
# 1/5e-3 = 200...
# 1/2e-3 = 500...

for dt = 1e-3
  @printf("Running with dt=%3.0e, fi=%.2f, ki=%d, and μ=%4.0e...\n", 
    dt, fi, ki, μ)

  plotname = @sprintf("steady_fi%02d_n%d_ki%02d_mu%4.0e_dt%4.0e", 
    Int(fi*10), n, ki, μ, dt)

  runforcingproblem(fi=fi, ki=ki, tf=tf, dt=dt, withplot=true, 
    ns=5, n=n, ν=ν, μ=μ, stepper="RK4", plotname=plotname,
    withoutput=true, filename=plotname, stochastic=false)
end
