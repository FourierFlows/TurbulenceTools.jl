module TwoDTurbTools

using TurbulenceTools, FourierFlows, FourierFlows.TwoDTurb, PyPlot, JLD2

import FourierFlows.TwoDTurb: energy, enstrophy, dissipation, work, drag, updatevars!

export cfl, getresidual, getdiags, savediags, restartchanproblem,
       runchanproblem, makechanproblem, initandrunchanproblem, makeplot,
       loadgridparams, loadtimestep, loadparams, loadlastsolution,
       loadforcingparams

"""
    loadgridparams(filename)

Returns nx, Lx, ny, Ly from the FourierFlows output stored in filename.
"""
function loadgridparams(filename)
  file = jldopen(filename)
  nx = file["grid/nx"]
  ny = file["grid/ny"]
  Lx = file["grid/Lx"]
  Ly = file["grid/Ly"]
  close(file)
  nx, Lx, ny, Ly
end

"""
    loadtimestep(filename)

Returns dt from the FourierFlows output stored in filename.
"""
function loadtimestep(filename)
  file = jldopen(filename)
  dt = file["timestepper/dt"]
  dt
end 

"""
    loadparams(filename)

Returns nu, nnu, mu, nmu from the FourierFlows output stored in filename.
"""
function loadparams(filename)
  file = jldopen(filename)
   nu = file["params/nu"]
  nnu = file["params/nnu"]
   mu = file["params/mu"]
  nmu = file["params/nmu"]
  close(file)
  nu, nnu, mu, nmu
end

"""
    loadforcingparams(filename)

Returns fi, ki from the FourierFlows output stored in filename.
"""
function loadforcingparams(filename)
  file = jldopen(filename)
  fi = file["forcingparams/fi"]
  ki = file["forcingparams/ki"]
  fi, ki
end

"""
    loadlastsolution(filename)

Returns the value of :sol with the highest timestep in the timeseris stored 
in the FourierFlows output file filename.
"""
function loadlastsolution(filename)
  file = jldopen(filename)
  laststep = parse(keys(file["timeseries/sol"])[end])
  sol = file["timeseries/sol/$laststep"]
  t = file["timeseries/t/$laststep"]
  laststep, t, sol
end

"""
    restartchanproblem(filename)

"""
function restartchanproblem(filename; ns=1, nt=1, tf=nothing, 
  stepper="FilteredRK4", withplot=false, plotname=nothing)

  fullfilename = filename[end-5:end] != ".jld2" ? filename*".jld2" : filename

  nx, Lx, ny, Ly = loadgridparams(fullfilename)
  nu, nnu, mu, nmu = loadparams(fullfilename)  
  fi, ki = loadforcingparams(fullfilename)
  dt = loadtimestep(fullfilename)
  step, t, sol = loadlastsolution(fullfilename) # Load last saved solution

  if tf != nothing
    nt = round(Int, (tf-t)/dt)
  else
    tf = t + nt*dt
  end

  prob, diags, nt = makechanproblem(n=nx, L=Lx, nu=nu, nnu=nnu, mu=mu, nmu=nmu, dt=dt, fi=fi, ki=ki, tf=tf, 
                                    stepper=stepper)
  prob.state.step = step
  prob.state.t = t
  prob.state.sol .= sol
  updatevars!(prob) # for fun

  newfilename = fullfilename[1:end-5] * "_restart"
  output = getbasicoutput(prob; filename=newfilename)
  runchanproblem(prob, diags, nt, fi; withplot=withplot, ns=ns, output=output, plotname=plotname)

  prob, diags, output
end

"""
    cfl(prob)

Returns the CFL number defined by CFL = max(u*dt/dX), where u = (U, V) is 
the horizontal velocity and dX = (dx, dy) the grid spacing.
"""
cfl(prob) = prob.ts.dt*maximum([maximum(prob.vars.U)/prob.grid.dx, maximum(prob.vars.V)/prob.grid.dy])

"""
    getresidual(prob, E, I, D, R, ψ, F; i0=1)

Returns the residual defined by

               dE
  residual  =  --  -  I  +  D  + R, 
               dt

where I = -<ψF>, D = nu<ψΔⁿζ>, and R = mu<ψΔⁿ¹ζ>, with n and n1 the order of 
hyper- and hypo-dissipation operators, respectively. For the stochastic case,
care is needed to calculate I correctly.
"""
function getresidual(prob, E, I, D, R, ψ, F, fi; ii0=1, iif=E.count) 

  # Forward difference: dEdt calculated at ii=ii0:(iif-1) 
  ii = ii0:(iif-1) 
  ii₊₁ = (ii0+1):iif

  # to calculate dEdt for fixed dt
  dEdt = ( E[ii₊₁] - E[ii] ) / prob.ts.dt
  dEdt - I[ii] + 0.5*fi^2 + D[ii] + R[ii]
end

function getresidual(prob, diags, fi; kwargs...)
  E, Z, D, I, R, F, ψ = diags[1:7]
  getresidual(prob, E, I, D, R, ψ, F, fi; kwargs...)
end

"""
    getdiags(prob, nt)

Returns a vector of Diagnostics that are meaningful in stochastically-forced
two-dimensional turbulence problems.
"""
function getdiags(prob, nt)
  forcing(prob) = deepcopy(prob.vars.Fh)

  getpsih(prob) = -prob.grid.invKKrsq.*prob.state.sol

  E = Diagnostic(energy,      prob, nsteps=nt)
  Z = Diagnostic(enstrophy,   prob, nsteps=nt)
  D = Diagnostic(dissipation, prob, nsteps=nt)
  I = Diagnostic(work,        prob, nsteps=nt)
  R = Diagnostic(drag,        prob, nsteps=nt)
  F = Diagnostic(forcing,     prob, nsteps=nt)
  ψ = Diagnostic(getpsih,     prob, nsteps=nt)

  [E, Z, D, I, R, F, ψ]
end

function savediags(out, diags)
  E, Z, D, I, R, F, ψ = diags[1:7]
  savediagnostic(E, "energy", out.filename)
  savediagnostic(Z, "enstrophy", out.filename)
  savediagnostic(D, "dissipation", out.filename)
  savediagnostic(I, "work", out.filename)
  savediagnostic(R, "drag", out.filename)
  savediagnostic(F, "forcing", out.filename)
  savediagnostic(ψ, "psih", out.filename)
  nothing
end

"""
    makeplot(prob, diags)

Makes a three-components plot of the vorticity field, energy tendency, and 
energy.
"""
function makeplot(prob, diags, fi; i₀=1)
  
  updatevars!(prob)  
  E, Z, D, I, R, F = diags

  close("all")
  fig, axs = subplots(ncols=3, nrows=1, figsize=(13, 4))

  sca(axs[1]); cla()
  pcolormesh(prob.grid.X, prob.grid.Y, prob.vars.q)
  xlabel(L"x")
  ylabel(L"y")

  sca(axs[2]); cla()
  plot(E[:time], -D[:], label="dissipation (\$D\$)")
  plot(E[:time], -R[:], label="drag (\$R\$)")
  
  ylabel("Energy sources and sinks")
  xlabel(L"t")
  legend(fontsize=10, loc="lower right")

  sca(axs[3]); cla()
  plot(E[:time], E[:])
  xlabel(L"t")
  ylabel(L"E")

  tight_layout()
  pause(0.01)
  nothing
end

"""
    initandrunchanproblem(; parameters...)

Create and run a two-dimensional turbulence problem with the "Chan forcing".
"""
function initandrunchanproblem(; n=128, L=2π, nu=4e-3, nnu=1, mu=1e-1, nmu=0, dt=1e-2, fi=1.0, ki=8, tf=10, ns=1, 
                               θ=π/4, withplot=false, withoutput=false, stepper="RK4", plotname=nothing,
                               filename="default")

  prob, diags, nt = makechanproblem(n=n, L=L, nu=nu, nnu=nnu, mu=mu, nmu=nmu, dt=dt, fi=fi, ki=ki, tf=tf, 
                                    stepper=stepper)
     
  if withoutput
    output = getbasicoutput(prob; filename=filename)
    jldopen(output.filename, "r+") do file
      file["forcingparams/fi"] = fi
      file["forcingparams/ki"] = ki
    end
  else
    output = nothing
  end

  runchanproblem(prob, diags, nt, fi; withplot=withplot, ns=ns, output=output, plotname=plotname)

  if withoutput; return prob, diags, output
  else;          return prob, diags
  end
end

"""
    prob, diags, nt = makechanproblem(; parameters...)

Returns a problem, vector of Diagnostics, and the number of timesteps for a 
two-dimensional turbulence problem forced by the "Chan forcing" 
with the specified parameters.
"""
function makechanproblem(; n=128, L=2π, nu=1e-3, nnu=1, mu=1e-1, nmu=-1, dt=1e-2, fi=1.0, ki=8, tf=1, stepper="RK4",
                         numdiags=Int(10^4))
  kii = ki*L/2π
  amplitude = fi*ki/sqrt(dt) * n^2/2
  function calcF!(F, sol, t, s, v, p, g)
    if t == s.t # not a substep
      F .= 0.0
      θk = 2π*rand() 
      phase = 2π*im*rand()
      i₁ = round(Int, abs(kii*cos(θk))) + 1
      j₁ = round(Int, abs(kii*sin(θk))) + 1  # j₁ >= 1
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
  prob = TwoDTurb.Problem(nx=n, Lx=L, nu=nu, nnu=nnu, mu=mu, nmu=nmu, dt=dt, calcF=calcF!, stepper=stepper)
  numdiags = minimum([nt, numdiags])
  diags = getdiags(prob, numdiags)

  prob, diags, nt
end

"""
    runchanproblem(prob, diags, nt, fi; ns=1, withplot=false, output=nothing, plotname=nothing, message=nothing)

Run the problem "prob" with useful messages, making plots if withplot 
and saving data if output is not nothing ns times. An optional additional 
message can be specified with the "message" keyword argument, where 
message(prob) is a function that returns a string. The string
"plotname" should be specified without the ".png" suffix.
"""
function runchanproblem(prob, diags, nt, fi; ns=1, withplot=false, output=nothing, plotname=nothing, message=nothing)

  @printf("\nRunning Chan 2D turbulence problem for %d steps...\n", nt)
  if output != nothing; println("Output: $(output.filename)"); end

  nint = round(Int, nt/ns)
  for i = 1:ns
    tc = @elapsed stepforward!(prob, diags, nint)
    updatevars!(prob)
    saveoutput(output)

    @printf("step: %04d, t: %.2e, cfl: %.3f, tc: %.2f s\n", prob.step, prob.t, cfl(prob), tc)

    if withplot     
      makeplot(prob, diags, fi)
      if plotname != nothing
        fullplotname = @sprintf("%s_%09d.png", plotname, prob.step)
        savefig(fullplotname, dpi=240)
      end
    end
  end

  updatevars!(prob)
  nothing
end

end # module
