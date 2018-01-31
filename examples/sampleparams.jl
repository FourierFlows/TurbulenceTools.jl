# Sample parameters

  n = 512
nki = 16

Ro = 0.2      # Rossby number
Rμ = 1e-4     # Ratio of drag and inertial time-scale
Rν = 1e-1     # Ratio of viscous and inertial time-scale
 f = 1e-4     # Planetary vorticity

tq = 1/(f*Ro) # Eddy turnover time-scale
tμ = 1/(f*Rμ) # Bottom drag time-scale
tν = 1/(f*Rν) # Viscous time-scale

 L = 1600e3
 μ = 1/tμ
kν = 2π/L * n/3 
nν = 2
 ν = 1/(kν^(2nν)*tν)
ki = nki*2π/L
fi = f*Ro/ki * sqrt(μ) # P = fi²

dt = tq/5
tf = 5tμ
nt = round(Int, tf/dt)
ns = 10
