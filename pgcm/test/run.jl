include("kernel.jl")

# path to folder containing param.h5 and <time step>.h5
path = string("runs/",ARGS[1])

# time step to start from
i0 = parse(Int, ARGS[2])

# output saving frequency (in timesteps)
save_freq = 10

# load model setup
println("loading setup")
m = load(path)

# load model state
println("loading state")
s = load(path, i0)

# time stepping
while s.i < m.i1
  # time step
  timestep!(m, s)
  # save every save_freq time steps
  if s.i%save_freq == 0
    save(s, path)
    println("saved")
  end
end
