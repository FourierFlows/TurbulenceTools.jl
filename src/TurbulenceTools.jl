__precompile__()

"""
Module structure:

  TurbulenceTools
    ... general utilities

    TwoDTurbTools
      ... twodturb utilities
      InitialValueProblems
      SteadyForcingProblems
      StochasticForcingProblems

    VerticallyCosineBoussinesqTools
    ...

"""
module TurbulenceTools
using FourierFlows, PyPlot

export makesquare!, ticksoff!, getbasicoutput

# Stuff to help with plotting
makesquare!(ax) = ax[:set_aspect](1, adjustable="box")
makesquare!(axs::AbstractArray) = for ax in axs; makesquare!(ax); end

ticksoff!(a) = a[:tick_params](bottom=false, left=false, labelbottom=false, 
  labelleft=false)
ticksoff!(axs::AbstractArray) = for ax in axs; ticksoff!(ax); end

"""
    getbasicoutput(prob, filename="default")

Returns Output whose only field is the solution.
"""
function getbasicoutput(prob; filename="default")
  filedir = joinpath(".", "data")
  if !isdir(filedir); mkdir(filedir); end
  filename = joinpath(filedir, filename)
  getsol(prob) = deepcopy(prob.state.sol)
  Output(prob, filename, (:sol, getsol))
end

# ----------------------------------------------------------------------------- 
# TwoDTurbTools --------------------------------------------------------------- 
# ----------------------------------------------------------------------------- 
module TwoDTurbTools
using FourierFlows.TwoDTurb
export cfl

"""
    cfl(prob)

Returns the CFL number defined by CFL = max([max(U)*dx/dt max(V)*dy/dt]).
"""
function cfl(prob)
  prob.ts.dt*maximum(
    [maximum(prob.vars.U)/prob.grid.dx, maximum(prob.vars.V)/prob.grid.dy])
end

include(joinpath("twodturb", "stochasticforcing.jl"))
include(joinpath("twodturb", "steadyforcing.jl"))

end # TwoDTurbTools module


# ----------------------------------------------------------------------------- 
# TwoDTurbTools --------------------------------------------------------------- 
# ----------------------------------------------------------------------------- 
module VerticallyCosineTools
using FourierFlows.VerticallyCosineBoussinesq
export cfl

"""
    cfl(prob)

Returns the CFL number defined by CFL = max([max(U)*dx/dt max(V)*dy/dt]).
"""
function cfl(prob)
  prob.ts.dt*maximum(
    [maximum(prob.vars.U)/prob.grid.dx, maximum(prob.vars.V)/prob.grid.dy,
     maximum(prob.vars.u)/prob.grid.dx, maximum(prob.vars.v)/prob.grid.dy  ])
end

include(joinpath("verticallycosineboussinesq", "stochasticwaveturb.jl"))

end

end # module
