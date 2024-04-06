
#=function setInput()::userinputs
  Nx = 2048 #2048
  Nt = 1024 # 1024
  cry = 4 #2->ZnTe; 4->GaAs
  sigma_t = 1e-12
  sigma_x = 1e-3
  lambda0 = 10.6e-6
  I0 = 60e13
  str_prefix = "/home/illesg/cst/2d-calculations/"
  STR = str_prefix * "gpu-test-GaAs-$(I0/1e13)GW"
  #gamma = deg2rad(22)
  gamma = acos(ngp(lambda0, 300, cry) / nTHzo(0.5e12 * 2 * pi, 300, cry))
  dz = 1e-6
  z_end = 8.0e-3 + dz
  x = range(-sigma_x, sigma_x, Nx) * 5
  t = range(-sigma_t, sigma_t, Nt) * 20
  return userinputs(Nx, Nt, cry, sigma_t, sigma_x, lambda0, I0, STR, gamma, dz, z_end, x, t)
end =#
@kwdef struct userinputs
  Nx::Int = 2048
  Nt::Int = 1024
  cry::Int = 2 # 2->ZnTe; 4->GaAs
  sigma_t::Float64 = 1.0e-12
  sigma_x::Float64 = 2.0e-3
  lambda0::Float64 = 10.6e-6
  I0::Float64 = 50e13
  STR::String = "/home/illesg/cst/2d-calculations/gpu-test-ZnTe-$(I0/1e13)GW-wide"
  gamma::Float64 = acos(ngp(lambda0, 300, cry) / nTHzo(0.5e12 * 2 * pi, 300, cry))
  dz::Float64 = 1e-6
  z_end::Float64 = 8e-3 + dz
  x::Vector = range(-sigma_x, sigma_x, Nx) * 5
  t::Vector = range(-sigma_t, sigma_t, Nt) * 20
end

