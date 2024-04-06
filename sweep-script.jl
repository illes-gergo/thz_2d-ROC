include("main.jl")
#=
function sweep()
  taus = [2.0, 1.75, 1.5, 1.25, 1, 0.75, 0.5] * 1e-12
  for tau in taus
    runcalc(tau)
    GC.gc()
  end
end
=#

function sweep()
  I0s = [90,100] * 1e13
  for I0 in I0s
    runcalc(I0)
  end
end
