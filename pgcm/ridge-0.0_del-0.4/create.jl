include("kernel.jl")

# experiment name (determines folder in runs/)
name = ARGS[1]

# aspect ratio (alpha)
a = .2

# friction coefficient (rho)
r = .1

# depth and its derivatives
d = .1
g(x) = 1-exp.(-x.^2 ./(2*d^2))
gp(x) = x/d^2 .*exp.(-x.^2 ./(2*d^2))
h(x,y) = g(x) .*g(1-x) .*g(1+y) .*g(1-y)
hx(x,y) = (gp(x) .*g(1-x)-g(x) .*gp(1-x)) .*g(1+y) .*g(1-y)
hy(x,y) = g(x) .*g(1-x) .*(gp(1+y) .*g(1-y)-g(1+y) .*gp(1-y))

# input file option
topofile = true

# diffusivity map
#   uniform diffusivity:
#     k(x,y,s) = 1.00e-2*ones(s)
#   bottom-intensified diffusivity:
#     k(x,y,s) = 1.00e-1*exp.(-(s+1) .*h(x,y)/.1)
#k(x,y,s) = 1.00e-1*exp.(-(s+1) .*h(x,y)/.1)
k(x,y,s) = 1.00e-2*ones(s)

# input k file option
kfile = true

# initial buoyancy profile
#   exponential stratification
#     del = .25
#     bi(z) = del*exp.(z/del)
#   constant stratification
#     bi(z) = z
#

del = .4
bi(z) = del*(exp.(z/del)-1)

#alpha = 0.0
#bi(z) = z + alpha*0.5*(z.^2+z)

# meridional restoring profile
#   no restoring:
#     c(y) = zeros(y)
#   restoring in "Southern Ocean":
#     c(y) = (1-tanh.((y+.5)/.1))/2
#   pick restoring constant 100x the average diffusivity
#c(y) = zeros(y)
c(y) = (1 - tanh.((y+.5)/.1))/2

# simulation length
T = 1

# number of grid points
nx = 100
ny = 200
ns = 20

# time step
dt = 2.0e-5
acc = 1.
dt = dt*acc

# set up model
m = ModelSetup(a, r, T, h, topofile, k, kfile, c, bi, nx, ny, ns, dt, acc)

# initialize model state
s = ModelState(m)

# save model setup and initial state
path = @sprintf("runs/%s", name)
mkpath(path)
save(m, path)
save(s, path)
