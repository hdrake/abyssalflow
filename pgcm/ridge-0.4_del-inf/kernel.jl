#Pkg.add("HDF5")
using HDF5
using Printf
using LinearAlgebra
using SparseArrays
using SuiteSparse
using Statistics

"""
    ModelSetup

Includes all information about the model setup, i.e. physical parameters, bathymetry
information, numerical parameters, fields needed in the calculation, and inversion matrices
for the implicit diffusion and the barotropic streamfunction calculation.
"""
struct ModelSetup
  # aspect ratio
  a::Float64

  # friction
  r::Float64

  # number of grid points
  nx::Int
  ny::Int
  ns::Int

  # time step
  dt::Float64
  acc::Float64

  # final time step number
  i1::Float64

  # grid spacing
  dx::Float64
  dy::Float64
  ds::Float64

  # coordinates at cell centers
  xc::Array{Float64,3}
  yc::Array{Float64,3}
  sc::Array{Float64,3}

  # coordinates at cell faces
  xf::Array{Float64,3}
  yf::Array{Float64,3}
  sf::Array{Float64,3}

  # diffusivity at various cell faces
  kfx::Array{Float64,3}
  kfy::Array{Float64,3}
  kfs::Array{Float64,3}

  # diffusivity at the surface (centers in c, y)
  k0c::Array{Float64,3}

  # restoring constant at cell centers
  cc::Array{Float64,3}

  # depth and its partial derivatives at cell centers
  hc::Array{Float64,3}
  hxc::Array{Float64,3}
  hyc::Array{Float64,3}

  # depth and its partial derivatives at x-faces
  hfx::Array{Float64,3}
  hxfx::Array{Float64,3}
  hyfx::Array{Float64,3}

  # depth and its partial derivatives at y-faces
  hfy::Array{Float64,3}
  hxfy::Array{Float64,3}
  hyfy::Array{Float64,3}

  # depth and its partial derivatives at cell nodes
  hn::Array{Float64,3}
  hxn::Array{Float64,3}
  hyn::Array{Float64,3}

  # initial buoyancy and restoring profile
  bic::Array{Float64,3}

  # array of sigma diffusion matrices
  A::Array{SuiteSparse.CHOLMOD.Factor{Float64}}

  # barotropic inversion matrix
  B::SuiteSparse.UMFPACK.UmfpackLU{Float64,Int64}

  function ModelSetup(a, r, T, h, topofile, k, kfile, c, bi, nx, ny, ns, dt, acc)

    # Constructs model with analytic functions giving depth h(x,y), diffusivity
    # k(x,y,s), restoring c(y), and initial buoyancy bi(z). All fields needed
    # for the calculation are deduced from these functions and passed to the
    # full constructor.

    # final time step
    i1 = convert(Int, round(T/dt*acc))

    # grid spacing (assuming domain 0<x<1, -1<y<1)
    dx = 1/nx
    dy = 2/ny
    ds = 1/ns
    
    # grid at cell centers
    xc = reshape(collect(dx/2:dx:1-dx/2), (1,1,nx))
    yc = reshape(collect(-1+dy/2:dy:1-dy/2), (1,ny,1))
    sc = reshape(collect(-1+ds/2:ds:-ds/2), (ns,1,1))
    
    # grid at cell faces
    xf = reshape(collect(dx:dx:1-dx), (1,1,nx-1))
    yf = reshape(collect(-1+dy:dy:1-dy), (1,ny-1,1))
    sf = reshape(collect(-1+ds:ds:-ds), (ns-1,1,1))
        
    if kfile == true
       # load k file
       file = h5open("k.h5", "r")

       # diffusivity at cell faces
       kfx = permutedims(repeat(read(file, "kfx"), outer = [1, 1, 1]),[3,2,1])
       kfy = permutedims(repeat(read(file, "kfy"), outer = [1, 1, 1]),[3,2,1])
       kfs = permutedims(repeat(read(file, "kfs"), outer = [1, 1, 1]),[3,2,1])

       # diffusivity at the surface (centers in x, y)
       k0c = permutedims(repeat(read(file, "k0c"), outer = [1, 1, 1]),[3,2,1])

    else
       # diffusivity at cell faces
       kfx = k.(xf,yc,sc)
       kfy = k.(xc,yf,sc)
       kfs = k.(xc,yc,sf)

       # diffusivity at the surface (centers in x, y)
       k0c = k.(xc,yc,zeros(1,ny,nx))
    end
    
    # restoring at cell centers
    cc = c.(yc)

    if topofile == true
       # load topography file
       file = h5open("topo.h5", "r")
       
       # depth and its partial derivatives at cell centers
       hc = permutedims(repeat(read(file, "hc"), outer = [1, 1, 1]),[3,2,1])
       hxc = permutedims(repeat(read(file, "hxc"), outer = [1, 1, 1]),[3,2,1])
       hyc = permutedims(repeat(read(file, "hyc"), outer = [1, 1, 1]),[3,2,1])
   
       # depth and its partial derivatives at x-faces
       hfx = permutedims(repeat(read(file, "hfx"), outer = [1, 1, 1]),[3,2,1])
       hxfx = permutedims(repeat(read(file, "hxfx"), outer = [1, 1, 1]),[3,2,1])
       hyfx = permutedims(repeat(read(file, "hyfx"), outer = [1, 1, 1]),[3,2,1])
   
       # depth and its partial derivatives at y-faces
       hfy = permutedims(repeat(read(file, "hfy"), outer = [1, 1, 1]),[3,2,1])
       hxfy = permutedims(repeat(read(file, "hxfy"), outer = [1, 1, 1]),[3,2,1])
       hyfy = permutedims(repeat(read(file, "hyfy"), outer = [1, 1, 1]),[3,2,1])
 
       # depth and its partial derivates at cell nodes
       hn = permutedims(repeat(read(file, "hn"), outer = [1, 1, 1]),[3,2,1])
       hxn = permutedims(repeat(read(file, "hxn"), outer = [1, 1, 1]),[3,2,1])
       hyn = permutedims(repeat(read(file, "hyn"), outer = [1, 1, 1]),[3,2,1])
    else
       # depth and its partial derivatives at cell centers
       hc = h(xc,yc)
       hxc = hx(xc,yc)
       hyc = hy(xc,yc)
       
       # depth and its partial derivatives at x-faces
       hfx = h(xf,yc)
       hxfx = hx(xf,yc)
       hyfx = hy(xf,yc)

       # depth and its partial derivatives at y-faces
       hfy = h(xc,yf)
       hxfy = hx(xc,yf)
       hyfy = hy(xc,yf)
       
       # depth and its partial derivates at cell nodes
       hn = h(xf,yf)
       hxn = hx(xf,yf)
       hyn = hy(xf,yf)
    end

    # initial buoyancy profile at cell centers
    bic = bi.(sc .*hc)

    # call full constructor to calculate inversion matrices
    new(a, r, nx, ny, ns, dt, acc, i1, dx, dy, ds, xc, yc, sc, xf, yf, sf, kfx, kfy,
        kfs, k0c, cc, hc, hxc, hyc, hfx, hxfx, hyfx, hfy, hxfy, hyfy, hn, hxn,
        hyn, bic)
  end
    
  function ModelSetup(a, r, nx, ny, ns, dt, acc, i1, dx, dy, ds, xc, yc, sc, xf, yf,
      sf, kfx, kfy, kfs, k0c, cc, hc, hxc, hyc, hfx, hxfx, hyfx, hfy, hxfy,
      hyfy, hn, hxn, hyn, bic) 

    # constructor that builds the inversion matrices
        
    # sigma diffusion matrices
    A = Array{SuiteSparse.CHOLMOD.Factor{Float64}}(undef, (ny, nx))
    K = kfs .*(1 .+a.^2 .*sf.^2 .*(hxc.^2 .+hyc.^2)) ./hc.^2
    K0 = k0c ./hc.^2
    for i = 1:nx
      for j = 1:ny
        o = -dt .*K[:,j,i] ./ds.^2
        d = ones(ns)
        d[1:ns-1] .+= dt .*K[:,j,i] ./ds.^2
        d[2:ns] .+= dt .*K[:,j,i] ./ds.^2
        d[ns] += 2 .*dt .*K0[1,j,i] ./ds.^2 # surface BC
        A[j,i] = factorize(spdiagm(1 => o, 0 => d, -1 => o))
      end
    end

    # barotropic inversion matrix
    idx(j,i) = (i-1)*(ny-1)+j
    B = spzeros((nx-1)*(ny-1),(nx-1)*(ny-1))
    for i = 1:nx-1
      for j = 1:ny-1

        # (r/h p_x)_x
        # not adjacent to land
        if hfy[1,j,i+1]!=0
          # nor boundary
          if i<nx-1
            B[idx(j,i),idx(j,i+1)] += r ./(hfy[1,j,i+1] .*dx.^2)
          end
          B[idx(j,i),idx(j,i)] -= r ./(hfy[1,j,i+1] .*dx.^2)
        end
        # not adjacent to land
        if hfy[1,j,i]!=0
          B[idx(j,i),idx(j,i)] -= r ./(hfy[1,j,i] .*dx.^2)
          # nor boundary
          if i>1
            B[idx(j,i),idx(j,i-1)] += r ./(hfy[1,j,i] .*dx.^2)
          end
        end

        # (r/h p_y)_y
        # not adjacent to land
        if hfx[1,j+1,i]!=0
          # nor boundary
          if j<ny-1
            B[idx(j,i),idx(j+1,i)] += r ./(hfx[1,j+1,i] .*dy.^2)
          end
          B[idx(j,i),idx(j,i)] -= r ./(hfx[1,j+1,i] .*dy.^2)
        end
        # not adjacent to land
        if hfx[1,j,i]!=0
          B[idx(j,i),idx(j,i)] -= r ./(hfx[1,j,i] .*dy.^2)
          # nor boundary
          if j > 1
            B[idx(j,i),idx(j-1,i)] += r ./(hfx[1,j,i] .*dy.^2)
          end
        end

        # -(y/h p_y)_x
        # not adjacent to boundary
        if i < nx-1 && j < ny-1
          # nor land
          if hn[1,j,i+1]!=0 
            B[idx(j,i),idx(j+1,i+1)] -= yf[j] ./(4 .*hn[1,j,i+1] .*dx .*dy)
          end
        end
        # not adjacent to boundary
        if i < nx-1 && j > 1
          # nor land
          if hn[1,j,i+1]!=0  
            B[idx(j,i),idx(j-1,i+1)] += yf[j] ./(4 .*hn[1,j,i+1] .*dx .*dy)
          end
        end
        # not adjacent to boundary
        if i > 1 && j < ny-1
          # nor land
          if hn[1,j,i-1]!=0
            B[idx(j,i),idx(j+1,i-1)] += yf[j] ./(4 .*hn[1,j,i-1] .*dx .*dy)
          end
        end
        # not adjacent to boundary
        if i > 1 && j > 1
          # nor land
          if hn[1,j,i-1]!=0
            B[idx(j,i),idx(j-1,i-1)] -= yf[j] ./(4 .*hn[1,j,i-1] .*dx .*dy)
          end
        end

        # (y/h p_x)_y
        # not adjacent to boundary 
        if i < nx-1 && j < ny-1
          # nor land
          if hn[1,j+1,i]!=0
            B[idx(j,i),idx(j+1,i+1)] += yf[j+1] ./(4 .*hn[1,j+1,i] .*dx .*dy)
          end
        end
        # not adjacent to boundary
        if i > 1 && j < ny-1
          # nor land
          if hn[1,j+1,i]!=0
            B[idx(j,i),idx(j+1,i-1)] -= yf[j+1] ./(4 .*hn[1,j+1,i] .*dx .*dy)
          end
        end
        # not adjacent to boundary
        if i < nx-1 && j > 1
          # nor land
          if hn[1,j-1,i]!=0
            B[idx(j,i),idx(j-1,i+1)] -= yf[j-1] ./(4 .*hn[1,j-1,i] .*dx .*dy)
          end
        end
        # not adjacent to boundary
        if i > 1 && j > 1
          # nor land
          if hn[1,j-1,i]!=0
            B[idx(j,i),idx(j-1,i-1)] += yf[j-1] ./(4 .*hn[1,j-1,i] .*dx .*dy)
          end
        end
      end
    end
    B = factorize(B)

    # Construct object.
    new(a, r, nx, ny, ns, dt, acc, i1, dx, dy, ds, xc, yc, sc, xf, yf, sf, kfx, kfy,
        kfs, k0c, cc, hc, hxc, hyc, hfx, hxfx, hyfx, hfy, hxfy, hyfy, hn, hxn,
	hyn, bic, A, B)
  end

end

mutable struct ModelState

  # This type includes all information to describe the model state, i.e. the
  # buoyancy field and the time step number. Other variables like velocities can
  # be determined from the buoyancy field.

  # buoyancy
  bc::Array{Float64,3}

  # time step number
  i::Int

end

# initialize with bc = bic
ModelState(m) = ModelState(m.bic, 0)

# derivatives at cell faces of field a given at cell centers
dxcf(m, a) = (a[:,:,2:m.nx] .-a[:,:,1:m.nx-1]) ./m.dx
dycf(m, a) = (a[:,2:m.ny,:] .-a[:,1:m.ny-1,:]) ./m.dy
dscf(m, a) = (a[2:m.ns,:,:] .-a[1:m.ns-1,:,:]) ./m.ds

# derivatives at cell centers of field a given at cell centers
function dxcc(m::ModelSetup, a)
  ax = zeros(size(a))
  ax[:,:,1] = (-a[:,:,3] .+4 .*a[:,:,2] .-3 .*a[:,:,1]) ./(2 .*m.dx)
  ax[:,:,2:m.nx-1] = (a[:,:,3:m.nx] .-a[:,:,1:m.nx-2]) ./(2 .*m.dx)
  ax[:,:,m.nx] = (3 .*a[:,:,m.nx] .-4 .*a[:,:,m.nx-1] .+a[:,:,m.nx-2]) ./(2 .*m.dx)
  return ax
end
function dycc(m::ModelSetup, a)
  ay = zeros(size(a))
  ay[:,1,:] = (-a[:,3,:] .+4 .*a[:,2,:] .-3*a[:,1,:]) ./(2 .*m.dy)
  ay[:,2:m.ny-1,:] = (a[:,3:m.ny,:] .-a[:,1:m.ny-2,:]) ./(2 .*m.dy)
  ay[:,m.ny,:] = (3 .*a[:,m.ny,:] .-4 .*a[:,m.ny-1,:] .+a[:,m.ny-2,:]) ./(2 .*m.dy)
  return ay
end
function dscc(m::ModelSetup, a, a0)
  as = zeros(size(a))
  as[1,:,:] = (-a[3,:,:] .+4*a[2,:,:] .-3*a[1,:,:]) ./(2 .*m.ds)
  as[2:m.ns-1,:,:] = (a[3:m.ns,:,:] .-a[1:m.ns-2,:,:]) ./(2 .*m.ds)
  as[m.ns,:,:] = (4 .*a0 .-3 .*a[m.ns,:,:] .-a[m.ns-1,:,:]) ./(3 .*m.ds)
  return as
end

# derivatives at cell centers of field a given at cell faces (2D and 3D)
function dxfc(m::ModelSetup, a::Array{Float64,3}, a0, a1)
  ax = zeros(m.ns,m.ny,m.nx)
  ax[:,:,1] = (a[:,:,1] .-a0) ./m.dx
  ax[:,:,2:m.nx-1] = (a[:,:,2:m.nx-1] .-a[:,:,1:m.nx-2]) ./m.dx
  ax[:,:,m.nx] = (a1 .-a[:,:,m.nx-1]) ./m.dx
  return ax
end
function dxfc(m::ModelSetup, a::Array{Float64,2}, a0, a1)
  my, mx = size(a)
  ax = zeros(my,m.nx)
  ax[:,1] = (a[:,1] .-a0) ./m.dx
  ax[:,2:m.nx-1] = (a[:,2:m.nx-1] .-a[:,1:m.nx-2]) ./m.dx
  ax[:,m.nx] = (a1 .-a[:,m.nx-1]) ./m.dx
  return ax
end
function dyfc(m::ModelSetup, a::Array{Float64,3}, a0, a1)
  ay = zeros(m.ns,m.ny,m.nx)
  ay[:,1,:] = (a[:,1,:] .-a0) ./m.dy
  ay[:,2:m.ny-1,:] = (a[:,2:m.ny-1,:] .-a[:,1:m.ny-2,:]) ./m.dy
  ay[:,m.ny,:] = (a1 .-a[:,m.ny-1,:]) ./m.dy
  return ay
end
function dyfc(m::ModelSetup, a::Array{Float64,2}, a0, a1)
  my, mx = size(a)
  ay = zeros(m.ny,mx)
  ay[1,:] = (a[1,:] .-a0) ./m.dy
  ay[2:m.ny-1,:] = (a[2:m.ny-1,:] .-a[1:m.ny-2,:]) ./m.dy
  ay[m.ny,:] = (a1 .-a[m.ny-1,:]) ./m.dy
  return ay
end

# interpolation from cell centers to cell faces
ixcf(m::ModelSetup, a) = (a[:,:,1:m.nx-1] .+a[:,:,2:m.nx]) ./2
iycf(m::ModelSetup, a) = (a[:,1:m.ny-1,:] .+a[:,2:m.ny,:]) ./2
iscf(m::ModelSetup, a) = (a[1:m.ns-1,:,:] .+a[2:m.ns,:,:]) ./2

function flux_div_step!(m::ModelSetup, bc, hFx, hFy, hFs)
  # Needs fluxes to be zero in hFs as well
  hFs[repeat(m.hc,outer = [m.ns-1, 1, 1]) .== 0.] .= 0.
  # 0/0 = NaN in Julia so need to overwrite hc so that fluxes vanish
  hc_finite = copy(m.hc)
  hc_finite[m.hc .== 0] .= 1.
  # flux divergence step
  bc[:,:,1:m.nx-1] .-= m.dt .*hFx ./(hc_finite[:,:,1:m.nx-1] .*m.dx)
  bc[:,:,2:m.nx] .+= m.dt .*hFx ./(hc_finite[:,:,2:m.nx] .*m.dx)
  bc[:,1:m.ny-1,:] .-= m.dt .*hFy ./(hc_finite[:,1:m.ny-1,:] .*m.dy)
  bc[:,2:m.ny,:] .+= m.dt .*hFy ./(hc_finite[:,2:m.ny,:] .*m.dy)
  bc[1:m.ns-1,:,:] .-= m.dt .*hFs ./(hc_finite .*m.ds)
  bc[2:m.ns,:,:] .+= m.dt .*hFs ./(hc_finite .*m.ds)
end

function diffusion!(m::ModelSetup, bc)
  # explicit interior fluxes (including cross terms in sigma)
  hFx = -m.a.^2 .*m.kfx .*(m.hfx .*dxcf(m,bc) .-m.sc .*m.hxfx .*dscc(m,ixcf(m,bc),0))
  hFy = -m.a.^2 .*m.kfy .*(m.hfy .*dycf(m,bc) .-m.sc .*m.hyfy .*dscc(m,iycf(m,bc),0))
  hFs = m.a.^2 .*m.kfs .*m.sf .*(m.hxc .*dxcc(m,iscf(m,bc))
      .+m.hyc .*dycc(m,iscf(m,bc)))
  flux_div_step!(m::ModelSetup, bc, hFx, hFy, hFs)
  # implicit sigma fluxes
  for i = 1:m.nx
    for j = 1:m.ny
      if m.hc[1,j,i] != 0
      	bc[:,j,i] .= m.A[j,i]\bc[:,j,i]
      end
    end
  end
end

function barotropic_flow(m::ModelSetup, bc)
  # RHS ("JEBAR" term)
  g = -sum(m.sc .*bc,dims=1) .*m.ds
  R = m.hxn .*ixcf(m,dycf(m,g)) .-m.hyn .*iycf(m,dxcf(m,g))
  # solve for streamfunction
  p = reshape(m.B\reshape(R,(m.nx-1)*(m.ny-1)),(m.ny-1,m.nx-1))
  # get velocities
  Ux = -dyfc(m,p,0,0) ./m.hfx[1,:,:]
  Uy = dxfc(m,p,0,0) ./m.hfy[1,:,:]
  return Ux, Uy
end

function thermal_wind(m::ModelSetup, bc)
  # Ux
  R = (-m.hfx .*iscf(m,m.yc .*dycc(m,ixcf(m,bc)) .+m.r .*dxcf(m,bc))
       .+m.sf .*(m.yc .*m.hyfx .+m.r .*m.hxfx) .*dscf(m,ixcf(m,bc)))
  Ux = zeros(m.ns,m.ny,m.nx-1)
  Ux[2:m.ns,:,:] = cumsum(R,dims=1)*m.ds ./(m.yc.^2 .+m.r^2)
  Ux  .-= mean(Ux,dims=1)
  # Uy
  R = (m.hfy .*iscf(m,m.yf .*dxcc(m,iycf(m,bc)) .-m.r .*dycf(m,bc))
       .+m.sf .*(-m.yf .*m.hxfy .+m.r .*m.hyfy) .*dscf(m,iycf(m,bc)))
  Uy = zeros(m.ns,m.ny-1,m.nx)
  Uy[2:m.ns,:,:] = cumsum(R,dims=1)*m.ds ./(m.yf.^2 .+m.r^2)
  Uy  .-= mean(Uy,dims=1)
  return Ux, Uy
end

function sigma_velocity(m::ModelSetup, Ux, Uy)
  # flux divergence of Ux, Uy
  R = -dxfc(m,m.hfx .*Ux,0,0)-dyfc(m,m.hfy .*Uy,0,0)
  # integrate up from bottom
  Us = cumsum(R[1:m.ns-1,:,:],dims=1)*m.ds ./m.hc
  return Us
end

function velocities(m::ModelSetup, bc)
  # get barotropic velocities
  Utx, Uty = barotropic_flow(m, bc)
  # get baroclinic velocities
  Ubx, Uby = thermal_wind(m, bc)
  # combine
  Ux = Ubx  .+ reshape(Utx, (1,m.ny,m.nx-1))
  Uy = Uby  .+ reshape(Uty, (1,m.ny-1,m.nx))
  # determine Us from continuity
  Us = sigma_velocity(m, Ux, Uy)
  return Ux, Uy, Us
end

function advection!(m::ModelSetup, bc)
  # get velocities
  Ux, Uy, Us = velocities(m, bc)
  # calculate advective fluxes
  hFx = m.hfx .*Ux .*ixcf(m,bc)
  hFy = m.hfy .*Uy .*iycf(m,bc)
  hFs = m.hc .*Us .*iscf(m,bc)
  # time step with flux divergence
  flux_div_step!(m, bc, hFx, hFy, hFs)
end

function restoring!(m::ModelSetup, bc)
  # restore b to bi
  bc[:] = m.bic  .+ (bc .-m.bic) .*exp.(-m.cc .*m.dt)
end

function timestep!(m::ModelSetup, s::ModelState)
  # diffusion step
  diffusion!(m, s.bc)
  # advection step
  advection!(m, s.bc)
  # restoring step
  restoring!(m, s.bc)
  # increase time step number
  s.i += convert(Int, m.acc)
  # print diagnostics
  if s.i%4000 == 0
    @printf("%6d %5.2e %10.5e\n", s.i, s.i .*m.dt ./m.acc, sum(s.bc .*m.hc) .*m.dx .*m.dy .*m.ds)
  end
end

function save(m::ModelSetup, path)
  # save model setup
  file = h5open(@sprintf("%s/param.h5", path), "w")
  write(file, "a", m.a)
  write(file, "r", m.r)
  write(file, "nx", m.nx)
  write(file, "ny", m.ny)
  write(file, "ns", m.ns)
  write(file, "dt", m.dt)
  write(file, "acc", m.acc)
  write(file, "i1", m.i1)
  write(file, "dx", m.dx)
  write(file, "dy", m.dy)
  write(file, "ds", m.ds)
  write(file, "xc", m.xc)
  write(file, "yc", m.yc)
  write(file, "sc", m.sc)
  write(file, "xf", m.xf)
  write(file, "yf", m.yf)
  write(file, "sf", m.sf)
  write(file, "kfx", m.kfx)
  write(file, "kfy", m.kfy)
  write(file, "kfs", m.kfs)
  write(file, "k0c", m.k0c)
  write(file, "cc", m.cc)
  write(file, "hc", m.hc)
  write(file, "hxc", m.hxc)
  write(file, "hyc", m.hyc)
  write(file, "hfx", m.hfx)
  write(file, "hxfx", m.hxfx)
  write(file, "hyfx", m.hyfx)
  write(file, "hfy", m.hfy)
  write(file, "hxfy", m.hxfy)
  write(file, "hyfy", m.hyfy)
  write(file, "hn", m.hn)
  write(file, "hxn", m.hxn)
  write(file, "hyn", m.hyn)
  write(file, "bic", m.bic)
  close(file)
end

function load(path)
  # load model setup
  file = h5open(@sprintf("%s/param.h5", path), "r")
  a = read(file, "a")
  r = read(file, "r")
  nx = read(file, "nx")
  ny = read(file, "ny")
  ns = read(file, "ns")
  dt = read(file, "dt")
  acc = read(file, "acc")
  i1 = read(file, "i1")
  dx = read(file, "dx")
  dy = read(file, "dy")
  ds = read(file, "ds")
  xc = read(file, "xc")
  yc = read(file, "yc")
  sc = read(file, "sc")
  xf = read(file, "xf")
  yf = read(file, "yf")
  sf = read(file, "sf")
  kfx = read(file, "kfx")
  kfy = read(file, "kfy")
  kfs = read(file, "kfs")
  k0c = read(file, "k0c")
  cc = read(file, "cc")
  hc = read(file, "hc")
  hxc = read(file, "hxc")
  hyc = read(file, "hyc")
  hfx = read(file, "hfx")
  hxfx = read(file, "hxfx")
  hyfx = read(file, "hyfx")
  hfy = read(file, "hfy")
  hxfy = read(file, "hxfy")
  hyfy = read(file, "hyfy")
  hn = read(file, "hn")
  hxn = read(file, "hxn")
  hyn = read(file, "hyn")
  bic = read(file, "bic")
  close(file)
  return ModelSetup(a, r, nx, ny, ns, dt, acc, i1, dx, dy, ds, xc, yc, sc, xf, yf,
      sf, kfx, kfy, kfs, k0c, cc, hc, hxc, hyc, hfx, hxfx, hyfx, hfy, hxfy,
      hyfy, hn, hxn, hyn, bic)
end

function save(s::ModelState, path)
  # save model state
  file = h5open(@sprintf("%s/%010d.h5", path, s.i), "w")
  write(file, "bc", s.bc)
  close(file)
end

function load(path, i)
# load model state
  file = h5open(@sprintf("%s/%010d.h5", path, i), "r")
  bc = read(file, "bc")
  close(file)
  return ModelState(bc, i)
end
