__precompile__()
module TurbulenceTools
using FourierFlows, PyPlot

# Stuff to help with plotting
makesquare!(ax) = ax[:set_aspect](1, adjustable="box")
makesquare!(axs::AbstractArray) = for ax in axs; makesquare!(ax); end

ticksoff!(a) = a[:tick_params](bottom=false, left=false, labelbottom=false, 
  labelleft=false)
ticksoff!(axs::AbstractArray) = for ax in axs; ticksoff!(ax); end



module TwoDTurbTools

using FourierFlows
import FourierFlows.TwoDTurb

include("twodturb/utils.jl")
include("twodturb/plotting.jl")
include("twodturb/problems.jl")

end # module

module VerticallyCosineTools


end # module

end # module
