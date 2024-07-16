using PlotlyJS, LazyGrids
include("fuggvenyek.jl")
include("valtozok.jl")

cry = 7

inputs = userinputs(cry=cry)


c0 = 3e8
d_eff = deff(inputs.cry)
e0 = 8.854187817e-12

SHG_SHIFT = floor(Int, inputs.Nt / 4)

tMax = inputs.t[end] - inputs.t[1]
dt = inputs.t[2] - inputs.t[1]

xMax = inputs.x[end] - inputs.x[1]
dx = inputs.x[2] - inputs.x[1]

dOmega = 2 * pi / tMax
omega_ = range(0, length=inputs.Nt, step=dOmega)

omegaMax = omega_[end] - omega_[1]
omega = omega_ .- omegaMax ./ 2

dkx = 2 * pi / xMax
kx = range(0, length=inputs.Nx, step=dkx) .- inputs.Nx / 2 * dkx
#kx0 = sin(inputs.gamma) ./ inputs.lambda0 .* 2 .* pi * neo(inputs.lambda0, 300, inputs.cry)

kxMax = kx[end] - kx[1]

omega_diff = round(Int, 2 * pi * c0 / inputs.lambda0 / dOmega)
omega0 = omega_diff .* dOmega
lambda0 = 2 * pi * c0 / omega0


(ct, cx) = ndgrid(inputs.t, inputs.x)
(comega_, ckx) = ndgrid(omega, kx)

comega = comega_ .+ omega0

comegaSHG = comega_ .+ 2 * omega0

comegaTHz = comega_ .- omega[1]


clambda = c0 ./ comega * 2 * pi

cLambdaSHG = c0 ./ comegaSHG * 2 * pi

n = neo(clambda, 300, inputs.cry)
ng = ngp(clambda, 300, inputs.cry)

khi_eff = 2 * deffTHz(inputs.cry, Omega=comegaTHz, nOmega0=neo(lambda0, 300, inputs.cry))

k_omega = n .* comega ./ c0
kx_omega = real.(k_omega .* sin(inputs.gamma))
kz_omega = real.(k_omega .* cos(inputs.gamma))

k_omegaSHG = neo(cLambdaSHG, 300, inputs.cry) .* comegaSHG ./ c0
kx_omegaSHG = real.(k_omegaSHG .* sin(inputs.gamma))
kz_omegaSHG = real.(k_omegaSHG .* cos(inputs.gamma))

nTHz = nTHzo(comegaTHz, 300, inputs.cry)

k_omegaTHz = nTHz .* comegaTHz ./ c0
kz_omegaTHz = real.(sqrt.(Complex.(k_omegaTHz .^ 2 - ckx .^ 2)))
l1 = Layout(title=attr(text="Pumpa fázis törésmutató"))
s1 = scatter(x=clambda[:,1], y=n[:,1])
l2 = Layout(title=attr(text="THz fázis törésmutató"))
s2 = scatter(x=comegaTHz[:,1]/2/pi, y=nTHz[:,1])
l3 = Layout(title=attr(text="Pumpa csoport törésmutató"))
s3 = scatter(x=clambda[:,1], y=ng[:,1])
println("gamma = $(rad2deg(inputs.gamma))")
display(plot(s1,l1))
display(plot(s2,l2))
display(plot(s3,l3))
