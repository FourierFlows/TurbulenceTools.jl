__precompile__()
module TurbulenceTools

module TwoDTurbTools

using FourierFlows
import FourierFlows.TwoDTurb

include("twodturb/utils.jl")
include("twodturb/plotting.jl")
include("twodturb/problems.jl")

end

end # module
