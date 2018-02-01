module TwoDTurbTools

using TurbulenceTools, FourierFlows, FourierFlows.TwoDTurb, PyPlot, JLD2

export cfl, getresidual, getdiags, savediags, restartchanproblem,
       runchanproblem, makechanproblem, initandrunchanproblem, makeplot

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

Returns ν, nν, μ, nμ from the FourierFlows output stored in filename.
"""
function loadparams(filename)
  file = jldopen(filename)
   ν = file["params/ν"]
  nν = file["params/nν"]
   μ = file["params/μ"]
  nμ = file["params/nμ"]
  close(file)
  ν, nν, μ, nμ
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
  ν, nν, μ, nμ = loadparams(fullfilename)  
  fi, ki = loadforcingparams(fullfilename)
  dt = loadtimestep(fullfilename)
  step, t, sol = loadlastsolution(fullfilename) # Load last saved solution

  if tf != nothing
    nt = round(Int, (tf-t)/dt)
  else
    tf = t + nt*dt
  end

  prob, diags, nt = makechanproblem(n=nx, L=Lx, ν=ν, nν=nν, μ=μ, nμ=nμ, dt=dt, 
    fi=fi, ki=ki, tf=tf, stepper=stepper)

  prob.state.step = step
  prob.state.t = t
  prob.state.sol .= sol
  updatevars!(prob) # for fun

  newfilename = fullfilename[1:end-5] * "_restart"
  output = getbasicoutput(prob; filename=newfilename)

  runchanproblem(prob, diags, nt; withplot=withplot, ns=ns, output=output,
    plotname=plotname)
  
  prob, diags, output
end

"""
    cfl(prob)

Returns the CFL number defined by CFL = max(u*dt/dX), where u = (U, V) is 
the horizontal velocity and dX = (dx, dy) the grid spacing.
"""
function cfl(prob)
  prob.ts.dt*maximum(
    [maximum(prob.vars.U)/prob.grid.dx, maximum(prob.vars.V)/prob.grid.dy])
end

"""
    getresidual(prob, E, I, D, R, ψ, F; i0=1)

Returns the residual defined by

               dE
  residual  =  --  -  I  +  D  + R, 
               dt

where I = -<ψF>, D = ν<ψΔⁿζ>, and R = μ<ψΔⁿ¹ζ>, with n and n1 the order of 
hyper- and hypo-dissipation operators, respectively. For the stochastic case,
care is needed to calculate I correctly.
"""
function getresidual(prob, E, I, D, R, ψ, F; ii0=1, iif=E.count) 

  # Forward difference: dEdt calculated at ii=ii0:(iif-1) 
  ii = ii0:(iif-1) 
  ii₊₁ = (ii0+1):iif

  # to calculate dEdt for fixed dt
  dEdt = ( E[ii₊₁] - E[ii] ) / prob.ts.dt
  dEdt - I[ii] + D[ii] + R[ii]
end

function getresidual(prob, diags; kwargs...)
  E, Z, D, I, R, F, ψ = diags[1:7]
  getresidual(prob, E, I, D, R, ψ, F; kwargs...)
end

"""
    getdiags(prob, nt)

Returns a vector of Diagnostics that are meaningful in stochastically-forced
two-dimensional turbulence problems.
"""
function getdiags(prob, nt)
  forcing(prob) = deepcopy(prob.vars.F)

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
function makeplot(prob, diags; i₀=1)
  
  updatevars!(prob)  
  E, Z, D, I, R, F = diags

  close("all")
  fig, axs = subplots(ncols=3, nrows=1, figsize=(13, 4))

  sca(axs[1]); cla()
  pcolormesh(prob.grid.X, prob.grid.Y, prob.vars.q)
  xlabel(L"x")
  ylabel(L"y")

  sca(axs[2]); cla()
  dEdt = (E[(i₀+1):E.count] - E[i₀:E.count-1])/prob.ts.dt
  ii = (i₀+1):E.count

  # dEdt = I - D - R?
  dEdt₁ = I[ii] - D[ii] - R[ii]
  residual = dEdt - dEdt₁

  plot(E.time[ii], dEdt,   label=L"E_t")
  plot(E.time[ii], -D[ii], label="dissipation (\$D\$)")
  plot(E.time[ii], -R[ii], label="drag (\$R\$)")
  #plot(E.time[ii], residual, "c^", markersize=0.5, label="residual")
  #plot(E.time[ii], I[ii], "o", markersize=0.5, label="injection (\$I\$)")

  ylabel("Energy sources and sinks")
  xlabel(L"t")
  legend(fontsize=10, loc="lower right")

  sca(axs[3]); cla()
  plot(E.time[ii], E[ii])
  xlabel(L"t")
  ylabel(L"E")

  tight_layout()
  nothing
end

"""
    runchanproblem(; parameters...)

Create and run a two-dimensional turbulence problem with the "Chan forcing".
"""
function initandrunchanproblem(; n=128, L=2π, ν=4e-3, nν=1, 
  μ=1e-1, nμ=0, dt=1e-2, fi=1.0, ki=8, tf=10, ns=1, θ=π/4, 
  withplot=false, withoutput=false, stepper="RK4", plotname=nothing,
  filename="default")

  prob, diags, nt = makechanproblem(n=n, L=L, ν=ν, nν=nν, μ=μ, nμ=nμ, dt=dt,
     fi=fi, ki=ki, tf=tf, stepper=stepper)
  
  if withoutput
    output = getbasicoutput(prob; filename=filename)
    jldopen(output.filename, "r+") do file
      file["forcingparams/fi"] = fi
      file["forcingparams/ki"] = ki
    end
  else
    output = nothing
  end

  runchanproblem(prob, diags, nt; withplot=withplot, ns=ns, output=output,
    plotname=plotname)

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
function makechanproblem(; n=128, L=2π, ν=1e-3, nν=1, 
  μ=1e-1, nμ=-1, dt=1e-2, fi=1.0, ki=8, tf=1, stepper="RK4")

  kii = ki*L/2π
  amplitude = fi*ki/sqrt(dt) * n^2/4
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
  prob = TwoDTurb.Problem(nx=n, Lx=L, ν=ν, nν=nν, μ=μ, nμ=nμ, dt=dt, 
                          calcF=calcF!, stepper=stepper)
  diags = getdiags(prob, nt)

  prob, diags, nt
end

"""
    runchanproblem(prob, diags, nt; ns=1, withplot=false, output=nothing,
                   plotname=nothing, message=nothing)

Run the problem "prob" with useful messages, making plots if withplot 
and saving data if output is not nothing ns times. An optional additional 
message can be specified with the "message" keyword argument, where 
message(prob) is a function that returns a string. The string
"plotname" should be specified without the ".png" suffix.
"""
function runchanproblem(prob, diags, nt; ns=1, withplot=false, output=nothing,
                        plotname=nothing, message=nothing)

  @printf("\nRunning Chan 2D turbulence problem for %d steps...\n", nt)
  if output != nothing; println("Output: $(output.filename)"); end

  nint = round(Int, nt/ns)
  for i = 1:ns
    tc = @elapsed stepforward!(prob, diags, nint)
    updatevars!(prob)  

    res = getresidual(prob, diags) # residual = dEdt - I + D + R

    # Some analysis
    E, Z, D, I, R, F = diags

    iavg = (length(res)-nint+1):length(res)

    avgI = mean(I[iavg])
    norm = maximum([ mean(abs.(D[iavg])), mean(abs.(R[iavg])) ])
    resnorm = mean(res[iavg])/norm

    @printf(
      "step: %04d, t: %.2e, cfl: %.3f, tc: %.2f s, <res>: %.3e, <I>: %.2f\n", 
      prob.step, prob.t, cfl(prob), tc, resnorm, avgI)

    if message != nothing; println(message(prob)); end

    if withplot     
      makeplot(prob, diags)
      if plotname != nothing
        fullplotname = @sprintf("%s_%09d.png", plotname, prob.step)
        savefig(fullplotname, dpi=240)
      end
    end

    if output != nothing
      saveoutput(output)
    end
  end

  updatevars!(prob)
  nothing
end

end # module
