__precompile__()
module TurbulenceTools

using FourierFlows, PyPlot
import FourierFlows.TwoDTurb
import FourierFlows.TwoDTurb: energy, enstrophy, dissipation, injection, drag

export getdiags, savediags, cfl, getsteadyforcingproblem, runwithmessage,
       getstochasticforcingproblem, makeplot, getresidual, getsimpleoutput,
       runsteadyforcingproblem, runstochasticforcingproblem


function cfl(prob)
  prob.ts.dt*maximum(
    [maximum(prob.vars.V)/prob.grid.dx, maximum(prob.vars.U)/prob.grid.dy])
end


function runstochasticforcingproblem(; n=128, L=2π, ν=4e-3, nν=1, 
  μ=1e-1, nμ=-1, dt=1e-2, fi=1.0, ki=8, tf=10, ns=1, withplot=false, 
  output=nothing, stepper="RK4", plotname=nothing)

  prob, diags, nt = getstochasticforcingproblem(n=n, L=L, ν=ν, nν=nν, μ=μ,
     nμ=nμ, dt=dt, fi=fi, ki=ki, tf=tf, stepper=stepper)

  if output != nothing
    out = getsimpleoutput(prob)
    runwithmessage(prob, diags, nt; withplot=withplot, ns=ns, output=out,
      plotname=plotname, forcing="stochastic")
  else
    runwithmessage(prob, diags, nt; withplot=withplot, ns=ns, 
      plotname=plotname, forcing="stochastic")
  end
  nothing
end


function runsteadyforcingproblem(; n=128, L=2π, ν=4e-3, nν=1, μ=1e-1, nμ=-1, 
  dt=1e-2, fi=1.0, ki=8, θ=π/4, tf=10, ns=1, withplot=false, output=nothing,
  stepper="RK4")

  prob, diags, nt = getsteadyforcingproblem(n=n, L=L, ν=ν, nν=nν, μ=μ, nμ=nμ,
    dt=dt, fi=fi, ki=ki, θ=θ, tf=tf, stepper=stepper)

  if output != nothing
    out = getsimpleoutput(prob)
    runwithmessage(prob, diags, nt; withplot=withplot, ns=ns, output=out,
      plotname=plotname)
  else
    runwithmessage(prob, diags, nt; withplot=withplot, ns=ns, 
      plotname=plotname)
  end
  nothing
end


function getresidual(prob, E, I, D, R; i₀=1)
  dEdt₀ = (E[(i₀+1):E.count] - E[i₀:E.count-1])/prob.ts.dt
  ii = (i₀+1):E.count
  dEdt₁ = I[ii] - D[ii] - R[ii] # dEdt = I - D - R
  dEdt₀ - dEdt₁
end

function getresidual(prob, diags; kwargs...)
  E, Z, D, I, R = diags
  getresidual(prob, E, I, D, R; kwargs...)
end


function getdiags(prob, nt)
  E = Diagnostic(energy,      prob, nsteps=nt)
  Z = Diagnostic(enstrophy,   prob, nsteps=nt)
  D = Diagnostic(dissipation, prob, nsteps=nt)
  R = Diagnostic(drag,        prob, nsteps=nt)
  I = Diagnostic(injection,   prob, nsteps=nt)
  diags = [E, Z, D, I, R]

  diags
end


function getstochasticforcingproblem(; n=128, L=2π, ν=1e-3, nν=1, 
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
  diags = getdiags(prob, nt)

  prob, diags, nt
end


function getchan2012prob(n, ν, ki; dt=1e-2, tf=1000)
  getstochasticforcingproblem(n=n, ν=ν, nν=1, μ=0, dt=dt, fi=1, ki=ki, tf=tf)
end 


function getsimpleoutput(prob)
  getsol(prob) = deepcopy(prob.state.sol)
  Output(prob, filename, (:sol, getsol))
end


function savediags(out, diags)
  E, Z, D, I, R = diags
  savediagnostic(E, "energy", out.filename)
  savediagnostic(Z, "enstrophy", out.filename)
  savediagnostic(D, "dissipation", out.filename)
  savediagnostic(I, "injection", out.filename)
  savediagnostic(R, "drag", out.filename)
  nothing
end


function getsteadyforcingproblem(; n=128, L=2π, ν=2e-3, nν=1, μ=1e-1, nμ=-1, 
  dt=1e-2, fi=1.0, ki=8, θ=π/4, tf=10, stepper="RK4")
  
  i₁ = round(Int, abs(ki*cos(θ))) + 1
  j₁ = round(Int, abs(ki*sin(θ))) + 1  # j₁ >= 1
  j₂ = n + 2 - j₁                       # e.g. j₁ = 1 => j₂ = nl+1

  amplitude = fi*ki * n^2/4
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

function makeplot(prob, diags; forcing="steady")

  TwoDTurb.updatevars!(prob)  
  E, Z, D, I, R = diags

  close("all")
  fig, axs = subplots(ncols=3, nrows=1, figsize=(13, 4))

  sca(axs[1]); cla()
  pcolormesh(prob.grid.X, prob.grid.Y, prob.vars.q)
  xlabel(L"x")
  ylabel(L"y")

  sca(axs[2]); cla()

  i₀ = 1
  dEdt = (E[(i₀+1):E.count] - E[i₀:E.count-1])/prob.ts.dt
  ii = (i₀+1):E.count

  # dEdt = I - D - R?
  total = I[ii] - D[ii] - R[ii]
  residual = dEdt - total

  plot(E.time[ii], I[ii], label="injection (\$I\$)")
  plot(E.time[ii], -D[ii], label="dissipation (\$D\$)")
  plot(E.time[ii], -R[ii], label="drag (\$R\$)")
  plot(E.time[ii], residual, "c-", label="residual")

  if forcing == "steady"
    plot(E.time[ii], total, label=L"I-D-R")
    plot(E.time[ii], dEdt, "k:", label=L"E_t")
  end

  ylabel("Energy sources and sinks")
  xlabel(L"t")
  legend(fontsize=10)

  sca(axs[3]); cla()
  plot(E.time[ii], E[ii])
  xlabel(L"t")
  ylabel(L"E")

  tight_layout()
  pause(0.1)

  nothing
end


function runwithmessage(prob, diags, nt; ns=1, withplot=false, output=nothing,
                        forcing="steady", plotname=nothing)
  for i = 1:ns
    tic()
    stepforward!(prob, diags, round(Int, nt/ns))
    tc = toq()
    TwoDTurb.updatevars!(prob)  
    res = getresidual(prob, diags) # residual = dEdt - I + D + R
    @printf("step: %04d, t: %.1f, cfl: %.3f, time: %.2f s, mean(res) = %.3e\n", 
      prob.step, prob.t, cfl(prob), tc, mean(res))

    if withplot     
      makeplot(prob, diags; forcing=forcing)
      if plotname != nothing
        fullplotname = @sprintf("%s_%d.png", plotname, prob.step)
        savefig(fullplotname, dpi=240)
      end
    end

    if output != nothing
      saveoutput(out)
    end

  end
end

end # module
