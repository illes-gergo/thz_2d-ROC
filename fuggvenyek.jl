using Symbolics

function deffTHz(cry; Omega::Array=[0], nOmega0::Float64=0)
  if cry == 0 # LN
    deff_ = 168e-12 .* ones(size(Omega))
  elseif cry == 2 # ZnTe
    dE = 4.2e-12
    C = -0.07
    omega0 = 2 * pi * 5.35e12
    Lambda0 = 2 * pi * 0.09e12
    r41 = @. dE * (1 + ((C * omega0^2) / (omega0^2 - Omega^2 - 1im * Lambda0 * Omega)))
    deff_ = @. abs(2 / sqrt(3) * (nOmega0^4 * r41) / (1 + C) / 4)
    deff_[isnan.(deff_)] .= 0
  elseif cry == 3 # GaP
    deff_ = 0
  elseif cry == 4 # GaAs
    deff_ = 2 / sqrt(3) * 86.5e-12 .* ones(size(Omega))
  elseif cry == 7 # ZnSe
    deff_ = 0
  end
  return deff_
end

function neo(lambda, T, cry)
  if cry == 4 #GaAs Skauli et al. 2003 0.97-17 um
    l = lambda * 1e6
    a0 = 4.372514
    a = [5.466742 0.0242996 1.957522]
    b = [0.4431307 0.8746453 36.9166]

    n = abs.(sqrt.(Complex.(a0 .+ 1 .+ a[1] * l .^ 2 ./ (l .^ 2 .- b[1]^2) + a[2] * l .^ 2 ./ (l .^ 2 .- b[2]^2) + a[(3)] * l .^ 2 ./ (l .^ 2 .- b[(3)]^2))))
    return n


  elseif cry == 0 #% LN
    lambda1 = lambda * 1e6
    a1 = 5.078
    a2 = 0.0964
    a3 = 0.2065
    a4 = 61.16
    a5 = 10.55
    a6 = 1.59e-2
    b1 = 4.677e-7
    b2 = 7.822e-8
    b3 = -2.653e-8
    b4 = 1.096e-4
    T = T - 273
    f = (T - 24.5) * (T + 570.82)

    n = @. real.(sqrt.(Complex.(a1 + b1 * f + (a2 + b2 * f) ./ (lambda1 .^ 2 - (a3 + b3 * f)^2) + (a4 + b4 * f) ./ (lambda1 .^ 2 - a5^2) - a6 * lambda1 .^ 2)))
    return n

  elseif cry == 2 # ZnTe
    n = @. real(sqrt(Complex(9.92 + 0.42530 / ((lambda * 1e6)^2 - 0.37766^2) + 8414.13 / ((lambda * 1e6)^2 - 56.5^2))))
  elseif cry == 7 # ZnSe
    a = 4.45813734
    b = -0.200859853
    c = 0.467216334
    d = -0.391371166
    e = 2.89566290
    f = -47.1362108
    n = @. real(sqrt(Complex(a * lambda1^2 / (lambda1^2 - b^2) + c * lambda1^2 / (lambda1^2 - d^2) + e * lambda1^2 / (lambda1^2 - f^2)) + 1))
    return n
  end
end

function ngp(lambda, T, cry)
  lambda1 = lambda * 1e6
  if cry == 4
    @variables l

    a0 = 4.372514
    a = [5.466742 0.0242996 1.957522]
    b = [0.4431307 0.8746453 36.9166]

    n0 = real.(sqrt.(a0 + 1 + a[(1)] * l .^ 2 ./ (l .^ 2 - b[(1)]^2) + a[(2)] * l .^ 2 ./ (l .^ 2 - b[(2)]^2) + a[(3)] * l .^ 2 ./ (l .^ 2 - b[(3)]^2)))

    a = n0 - l * Symbolics.derivative(n0, l)
    ng = Symbolics.value(substitute(a, l => lambda1))
    return ng

  elseif cry == 2
    @variables l
    n0 = real(sqrt((9.92 + 0.42530 / (l^2 - 0.37766^2) + 8414.13 / (l^2 - 56.5^2))))
    a = n0 - l * Symbolics.derivative(n0, l)
    ng = Symbolics.value(substitute(a, l => lambda1))
    return ng
  elseif cry == 0
    @variables l
    a1 = 5.078
    a2 = 0.0964
    a3 = 0.2065
    a4 = 61.16
    a5 = 10.55
    a6 = 1.59e-2
    b1 = 4.677e-7
    b2 = 7.822e-8
    b3 = -2.653e-8
    b4 = 1.096e-4
    T = T - 273
    f = (T - 24.5) * (T + 570.82)

    n0 = sqrt((a1 + b1 * f + (a2 + b2 * f) ./ (l .^ 2 - (a3 + b3 * f)^2) + (a4 + b4 * f) ./ (l .^ 2 - a5^2) - a6 * l .^ 2))
    a = n0 - l * Symbolics.derivative(n0, l)
    ng = Symbolics.value(substitute(a, l => lambda1))
    return ng
  end
end

function nTHzo(omega, T, cry)
  if cry == 0

    if T == 100
      nTHz = 4.69732 - 0.03006e-12 * omega / 2 / pi + 0.04066e-24 * (omega / 2 / pi) .^ 2
    elseif T == 300
      nTHz = @. 4.91372 - 0.01747e-12 * omega / 2 / pi + 0.04004e-24 * (omega / 2 / pi) .^ 2
    end


  elseif cry == 4 || cry == 2
    nTHz = real.(sqrt.(er(omega, T, cry)))
  end
  return nTHz
end

function n2value(cry)
  if cry == 4 # GaAs
    n2_ = 5.9e-18
  elseif cry == 2
    n2_ = 2.25e-18
  elseif cry == 0
    n2_ = 9.33e-20
  end
  return n2_
end

function deff(cry)
  if cry == 4 #% GaAs
    deff_ = 2 / sqrt(3) * 80e-12
  elseif cry == 2
    deff_ = 92.2e-12
  elseif cry == 0
    deff_ = 0
  end
  return deff_
end

function aTHzo(omega, T, cry)
  if cry != 0
    alpha = abs.(-2 .* omega / 3e8 .* imag(sqrt.(er(omega, T, cry))))
  elseif cry == 0
    alpha = @. 100 * (2.16411e-12 * omega / 2 / pi + 10.81e-24 * (omega / 2 / pi) .^ 2)
  end
  alpha[alpha.>1e5] .= 1e5
  return alpha
end

function er(omega, T, cry)
  if cry == 4 #GaAs
    if T == 300 #ord
      nu = omega / 2 / pi / 3e8 * 0.01
      e_inf = 11
      nu_L = 292.1
      nu_T = 268.7
      G = 2.4
      nu = omega / 2 / pi / 3e8 * 1e-2

      er_ = e_inf * (1 .+ (nu_L^2 .- nu_T^2) ./ (nu_T^2 .- nu .^ 2 .+ 1im * G * nu))
    end
  elseif cry == 2
    A = [4.262e-2, 1.193e-2, 3, 0.008, 0.0029, 0.005]
    oi = [56, 94.82, 182, 239, 299, 351]
    gi = [13.06, 8.18, 0.9, 6.7, 9.3, 10.3]
    es = 6.8

    tsc = 200e-15
    meff = 0.11 * 9.109e-31
    q = 1.602e-19
    e0 = 8.8541878e-12
    Nc = 0
    op2 = q^2 * Nc / e0 / meff

    er_ = @. es + A[1] * oi[1] .^ 2 ./ (oi[1] .^ 2 - omega .^ 2 * 0.01^2 / (2 * pi * 3e8)^2 - 2 * 1im * gi[1] * omega * 0.01 / (2 * pi * 3e8)) +
             A[2] * oi[2] .^ 2 ./ (oi[2] .^ 2 - omega .^ 2 * 0.01^2 / (2 * pi * 3e8)^2 - 2 * 1im * gi[2] * omega * 0.01 / (2 * pi * 3e8)) +
             A[3] * oi[3] .^ 2 ./ (oi[3] .^ 2 - omega .^ 2 * 0.01^2 / (2 * pi * 3e8)^2 - 2 * 1im * gi[3] * omega * 0.01 / (2 * pi * 3e8)) +
             A[4] * oi[4] .^ 2 ./ (oi[4] .^ 2 - omega .^ 2 * 0.01^2 / (2 * pi * 3e8)^2 - 2 * 1im * gi[4] * omega * 0.01 / (2 * pi * 3e8)) +
             A[5] * oi[5] .^ 2 ./ (oi[5] .^ 2 - omega .^ 2 * 0.01^2 / (2 * pi * 3e8)^2 - 2 * 1im * gi[5] * omega * 0.01 / (2 * pi * 3e8)) +
             A[6] * oi[6] .^ 2 ./ (oi[6] .^ 2 - omega .^ 2 * 0.01^2 / (2 * pi * 3e8)^2 - 2 * 1im * gi[6] * omega * 0.01 / (2 * pi * 3e8)) -
             op2 ./ (omega .^ 2 + 1im * omega / tsc)
    if typeof(omega) != Float64
      er_[isnan.(er_)] .= 1
      er_[isinf.(er_)] .= 1
    end
  end
  return er_
end
