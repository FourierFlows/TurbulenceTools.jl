__precompile__()
module TurbulenceTools
using FourierFlows
import FourierFlows.TwoDTurb

include("utils.jl")
include("plotting.jl")
include("problems.jl")

end # module
