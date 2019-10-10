include("kernel.jl")

# path to folder containing param.h5 and <time step>.h5
path = string("runs/",args[1])

# time step to start from
i0 = parse(Int, args[2])

# output saving frequency (in timesteps)
save_freq = 5000

# load model setup
m = load(path)

# load model state
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
