export makesteadyforcingproblem, runwithmessage, makestochasticforcingproblem,
       runforcingproblem

"""
    runforcingproblem(; parameters...)

Create and run a two-dimensional turbulence problem with the "Chan forcing".
"""
function runforcingproblem(; n=128, L=2π, ν=4e-3, nν=1, 
  μ=1e-1, nμ=-1, dt=1e-2, fi=1.0, ki=8, tf=10, ns=1, θ=π/4, 
  withplot=false, withoutput=false, stepper="RK4", plotname=nothing,
  filename="default", stochastic=false)

  if stochastic
    prob, diags, nt = makestochasticforcingproblem(n=n, L=L, ν=ν, nν=nν, μ=μ,
       nμ=nμ, dt=dt, fi=fi, ki=ki, tf=tf, stepper=stepper)
  else
    prob, diags, nt = makesteadyforcingproblem(n=n, L=L, ν=ν, nν=nν, μ=μ, nμ=nμ,
      dt=dt, fi=fi, ki=ki, θ=θ, tf=tf, stepper=stepper)
  end

  if withoutput
    output = getbasicoutput(prob; filename=filename)
    runwithmessage(prob, diags, nt; withplot=withplot, ns=ns, output=output,
      plotname=plotname, stochasticforcing=stochastic)
  else
    runwithmessage(prob, diags, nt; withplot=withplot, ns=ns, 
      plotname=plotname, stochasticforcing=stochastic)
  end

  nothing
end


function makestochasticforcingproblem(; n=128, L=2π, ν=1e-3, nν=1, 
  μ=1e-1, nμ=-1, dt=1e-2, fi=1.0, ki=8, tf=1, stepper="RK4")

  amplitude = fi*ki/sqrt(dt) * n^2/4
  function calcF!(F, sol, t, s, v, p, g)
    if t == s.t # not a substep
      F .= 0.0
      θk = 2π*rand() 
      phase = 2π*im*rand()
      i₁ = round(Int, abs(ki*cos(θk))) + 1
      j₁ = round(Int, abs(ki*sin(θk))) + 1  # j₁ >= 1
      j₂ = g.nl + 2 - j₁                    # e.g. j₁ = 1 => j₂ = nl+1
      if j₁ != 1  # apply forcing to l = (+/-)l★ mode
        F[i₁, j₁] = amplitude*exp(phase)
        F[i₁, j₂] = amplitude*exp(phase)
      else        # apply forcing to l=0 mode
        F[i₁, 1] = 2amplitude*exp(phase)
      end
    end
    nothing
  end

  nt = round(Int, tf/dt)
  prob = TwoDTurb.ForcedProblem(nx=n, Lx=L, ν=ν, nν=nν, μ=μ, nμ=nμ, dt=dt, 
    calcF=calcF!, stepper=stepper)
  diags = getdiags(prob, nt; stochasticforcing=true)

  prob, diags, nt
end


function getchan2012prob(n, ν, ki; dt=1e-2, tf=1000)
  makestochasticforcingproblem(n=n, ν=ν, nν=1, μ=0, dt=dt, fi=1, ki=ki, tf=tf)
end 


function makesteadyforcingproblem(; n=128, L=2π, ν=2e-3, nν=1, μ=1e-1, nμ=-1, 
  dt=1e-2, fi=1.0, ki=8, θ=π/4, tf=10, stepper="RK4")
  
  i₁ = round(Int, abs(ki*cos(θ))) + 1
  j₁ = round(Int, abs(ki*sin(θ))) + 1  # j₁ >= 1
  j₂ = n + 2 - j₁                      # e.g. j₁ = 1 => j₂ = nl+1
  amplitude = fi*ki * n^2/4

  # F = fi*ki*cos(i)*cos(j) (essentially)
  function calcF!(F, sol, t, s, v, p, g)
    if s.step == 1
      F[i₁, j₁] = amplitude
      F[i₁, j₂] = amplitude
    end
    nothing
  end

  nt = round(Int, tf/dt)
  prob = TwoDTurb.ForcedProblem(nx=n, Lx=L, ν=ν, nν=nν, μ=μ, nμ=nμ, dt=dt, 
    calcF=calcF!, stepper=stepper)
  TwoDTurb.set_q!(prob, rand(prob.grid.nx, prob.grid.ny))
  diags = getdiags(prob, nt)

  prob, diags, nt
end


function runwithmessage(prob, diags, nt; ns=1, withplot=false, output=nothing,
                        stochasticforcing=false, plotname=nothing)

  nint = round(Int, nt/ns)
  for i = 1:ns
    tic()
    stepforward!(prob, diags, nint)
    tc = toq()
    TwoDTurb.updatevars!(prob)  

    res = getresidual(prob, diags) # residual = dEdt - I + D + R

    # Some analysis
    E, Z, D, I, R, F = diags

    iavg = (length(res)-nint+1):length(res)

    avgI = mean(I[iavg])
    norm = maximum([ mean(abs.(D[iavg])), mean(abs.(R[iavg])) ])
    resnorm = mean(res[iavg])/norm

    @printf(
      "step: %04d, t: %.1f, cfl: %.3f, tc: %.2f s, <res>: %.3e, <I>: %.2f\n", 
      prob.step, prob.t, cfl(prob), tc, resnorm, avgI)

    if withplot     
      makeplot(prob, diags; stochasticforcing=stochasticforcing)
      if plotname != nothing
        plotdir = joinpath(".", "plots")
        fullplotname = joinpath(plotdir, 
          @sprintf("%s_%d.png", plotname, prob.step))
        if !isdir(plotdir); mkdir(plotdir); end
        savefig(fullplotname, dpi=240)
      end
    end

    if output != nothing
      saveoutput(output)
    end

  end
end

