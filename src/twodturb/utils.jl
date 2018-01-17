import FourierFlows.TwoDTurb: energy, enstrophy, dissipation, injection, drag

export cfl, getresidual, getdiags, getbasicoutput, savediags

"""
    cfl(prob)

Returns the CFL number defined by CFL = max([max(U)*dx/dt max(V)*dy/dt]).
"""
function cfl(prob)
  prob.ts.dt*maximum(
    [maximum(prob.vars.V)/prob.grid.dx, maximum(prob.vars.U)/prob.grid.dy])
end


"""
    getresidual(prob, E, I, D, R, F; i0=1)

Returns the residual defined by

               dE
  residual  =  --  -  I  +  D  + R, 
               dt

where I = -<ψF>, D = ν<ψΔⁿζ>, and R = μ<ψΔⁿ¹ζ>, with n and n1 the order of 
hyper- and hypo-dissipation operators, respectively. For the stochastic case,
care is needed to calculate the dissipation correctly.
"""
function getresidual(prob, E, I, D, R, ψ, F; ii0=1, iif=E.count, 
  stochasticforcing=false)

  # Forward difference: dEdt calculated at ii=ii0:(iif-1) 
  ii = ii0:(iif-1) 
  iifwd = (ii0+1):iif

  # to calculate dEdt for fixed dt
  dEdt = ( E[iifwd] - E[ii] ) / prob.ts.dt

  if stochasticforcing
    injkernel = 0.5*ψ[ii]*conj(F[ii]+F[iifwd])
    inj = -1/(g.Lx*g.Ly)*FourierFlows(injkernel)
  else
    inj = I[ii]
  end
  
  dEdt - inj + D[ii] + R[ii]
end

function getresidual(prob, diags; kwargs...)
  E, Z, D, I, R, F, ψ = diags[1:7]
  getresidual(prob, E, I, D, R, F; kwargs...)
end


function getdiags(prob, nt; stochasticforcing=false)
  forcing(prob) = deepcopy(prob.vars.F)

  getpsih(prob) = -prob.grid.invKKrsq.*prob.state.sol

  E = Diagnostic(energy,      prob, nsteps=nt)
  Z = Diagnostic(enstrophy,   prob, nsteps=nt)
  D = Diagnostic(dissipation, prob, nsteps=nt)
  I = Diagnostic(injection, prob, nsteps=nt)
  R = Diagnostic(drag,        prob, nsteps=nt)
  F = Diagnostic(forcing,     prob, nsteps=nt)
  ψ = Diagnostic(getpsi,      prob, nsteps=nt)

  [E, Z, D, I, R, F, ψ]
end


function getbasicoutput(prob; filename="default")
  filedir = joinpath(".", "data")
  if !isdir(filedir); mkdir(filedir); end
  filename = joinpath(filedir, filename)
  getsol(prob) = deepcopy(prob.state.sol)
  Output(prob, filename, (:sol, getsol))
end


function savediags(out, diags)
  E, Z, D, I, R, F, ψ = diags[1:7]
  savediagnostic(E, "energy", out.filename)
  savediagnostic(Z, "enstrophy", out.filename)
  savediagnostic(D, "dissipation", out.filename)
  savediagnostic(I, "injection", out.filename)
  savediagnostic(R, "drag", out.filename)
  savediagnostic(F, "forcing", out.filename)
  savediagnostic(ψ, "psih", out.filename)
  nothing
end
