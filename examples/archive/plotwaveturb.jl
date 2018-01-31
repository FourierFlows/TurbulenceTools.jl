using PyPlot, JLD2, FourierFlows, FourierFlows.VerticallyCosineBoussinesq

n = 256

filename = "testwaveturb_1.jld2"
file = jldopen(filename, "r")

steps = parse.(Int, keys(file["timeseries/sol"]))
m = file["params/m"]
N = file["params/N"]
nx = file["grid/nx"]
ny = file["grid/ny"]
Lx = file["grid/Lx"]
Ly = file["grid/Ly"]
close(file)

g = TwoDGrid(nx, Lx, ny, Ly)

fig, ax = subplots()
for step in steps

  try
    file = jldopen(filename, "r")
    sol = file["timeseries/sol/$step"]
    #@views Z = irfft(sol[:, :, 1], n)
    @views Zh = sol[:, :, 1]
    @views uh = sol[:, :, 2]
    @views vh = sol[:, :, 3]
    @views ph = sol[:, :, 4]
    Q = mode0apv(uh, vh, ph, Zh, m, N, g)

    cla(); imshow(Q); pause(0.1)
  finally
    close(file)
  end

end
