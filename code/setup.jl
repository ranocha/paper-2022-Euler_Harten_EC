
using Pkg; Pkg.activate(@__DIR__); Pkg.instantiate()

using Trixi

###############################################################################
# Implementation of the flux of Sjögreen and Yee (2019) for the compressible
# Euler equations on the 2D `TreeMesh` of Trixi.jl.
# https://doi.org/10.1007/s10915-019-01013-1

struct HartensCompressibleEulerEquations2D{RealT<:Real} <: Trixi.AbstractCompressibleEulerEquations{2, 4}
  equations::CompressibleEulerEquations2D{RealT}
  alpha::RealT # parameter of Harten's entropies
end

function HartensCompressibleEulerEquations2D(gamma::Real, alpha::Real)
  equations = CompressibleEulerEquations2D(gamma)
  HartensCompressibleEulerEquations2D(equations, alpha)
end


Trixi.varnames(f, equations::HartensCompressibleEulerEquations2D) = Trixi.varnames(f, equations.equations)

Trixi.flux(u, orientation_or_normal_direction, equations::HartensCompressibleEulerEquations2D) = flux(u, orientation_or_normal_direction, equations.equations)


Trixi.flux_shima_etal(u_ll, u_rr, orientation_or_normal_direction, equations::HartensCompressibleEulerEquations2D) = flux_shima_etal(u_ll, u_rr, orientation_or_normal_direction, equations.equations)

Trixi.flux_ranocha(u_ll, u_rr, orientation_or_normal_direction, equations::HartensCompressibleEulerEquations2D) = flux_ranocha(u_ll, u_rr, orientation_or_normal_direction, equations.equations)


Trixi.initial_condition_density_wave(x, t, equations::HartensCompressibleEulerEquations2D) = initial_condition_density_wave(x, t, equations.equations)

Trixi.initial_condition_weak_blast_wave(x, t, equations::HartensCompressibleEulerEquations2D) = initial_condition_weak_blast_wave(x, t, equations.equations)

Trixi.initial_condition_convergence_test(x, t, equations::HartensCompressibleEulerEquations2D) = initial_condition_convergence_test(x, t, equations.equations)

Trixi.source_terms_convergence_test(u, x, t, equations::HartensCompressibleEulerEquations2D) = source_terms_convergence_test(u, x, t, equations.equations)


Trixi.cons2prim(u, equations::HartensCompressibleEulerEquations2D) = cons2prim(u, equations.equations)

Trixi.max_abs_speeds(u, equations::HartensCompressibleEulerEquations2D) = Trixi.max_abs_speeds(u, equations.equations)

Trixi.max_abs_speed_naive(u_ll, u_rr, orientation_or_normal_direction, equations::HartensCompressibleEulerEquations2D) = Trixi.max_abs_speed_naive(u_ll, u_rr, orientation_or_normal_direction, equations.equations)


function Trixi.entropy(u, equations::HartensCompressibleEulerEquations2D)
  alpha = equations.alpha
  gamma = equations.equations.gamma
  rho, v1, v2, p = cons2prim(u, equations)

  return -(gamma + alpha) / (gamma - 1) * rho * (p / rho^gamma)^(1 / (alpha + gamma))
end

function Trixi.cons2entropy(u, equations::HartensCompressibleEulerEquations2D)
  alpha = equations.alpha
  gamma = equations.equations.gamma
  rho, v1, v2, p = cons2prim(u, equations)

  h = (p / rho^gamma)^(1 / (alpha + gamma))
  factor = rho / p * h
  w1 = -alpha / (gamma - 1) * h - 0.5 * factor * (v1^2 + v2^2)
  w2 = factor * v1
  w3 = factor * v2
  w5 = -factor

  return SVector(w1, w2, w3, w5)
end

function flux_potential(u, orientation::Integer, equations::HartensCompressibleEulerEquations2D)
  alpha = equations.alpha
  gamma = equations.equations.gamma
  rho, v1, v2, p = cons2prim(u, equations)

  h = (p / rho^gamma)^(1 / (alpha + gamma))

  if orientation == 1
    return rho * h * v1
  else # if orientation == 2
    return rho * h * v2
  end
end


function cons2param(u, equations::HartensCompressibleEulerEquations2D)
  alpha = equations.alpha
  gamma = equations.equations.gamma
  rho, v1, v2, p = cons2prim(u, equations)

  h = (p / rho^gamma)^(1 / (alpha + gamma))
  factor = rho / p * h
  z1 = factor
  z2 = v1
  z3 = v2
  z5 = p

  return SVector(z1, z2, z3, z5)
end

function exponential_average(a, b, beta_minus_1)
  beta = beta_minus_1 + 1

  r = a / b
  if abs(r - 1) < 1.0e-4
    return b^beta_minus_1 * (
      1 + (beta - 1) * (r - 1) / 2 + (beta - 1) * (beta - 2) * (r - 1)^2 / 6
    )
  else
    return (a^beta - b^beta) / (beta * (a - b))
  end
end

function flux_sjögreen_yee(u_ll, u_rr, orientation::Integer, equations::HartensCompressibleEulerEquations2D)
  alpha = equations.alpha
  gamma = equations.equations.gamma

  # Unpack left and right state
  z1_ll, z2_ll, z3_ll, z5_ll = cons2param(u_ll, equations)
  z1_rr, z2_rr, z3_rr, z5_rr = cons2param(u_rr, equations)

  if orientation == 1
    zhat_ll = z2_ll
    zhat_rr = z2_rr
  else # if orientation == 2
    zhat_ll = z3_ll
    zhat_rr = z3_rr
  end

  z1_zhat_mean = 0.5 * (z1_ll * zhat_ll + z1_rr * zhat_rr)
  z1_power_mean = 0.5 * (z1_ll^(-gamma / alpha) + z1_rr^(-gamma / alpha))
  z5_power_expmean = exponential_average(z5_ll, z5_rr, (1 - gamma - alpha) / alpha)
  h1 = z1_zhat_mean / (z1_power_mean * z5_power_expmean)

  if orientation == 1
    h2 = h1 * 0.5 * (z2_ll + z2_rr) + 0.5 * (z5_ll + z5_rr)
    h3 = h1 * 0.5 * (z3_ll + z3_rr)
  else # if orientation == 2
    h2 = h1 * 0.5 * (z2_ll + z2_rr)
    h3 = h1 * 0.5 * (z3_ll + z3_rr) + 0.5 * (z5_ll + z5_rr)
  end

  z5_power_mean = 0.5 * (z5_ll^(-(gamma - 1) / alpha) + z5_rr^(-(gamma - 1) / alpha))
  z1_power_expmean = exponential_average(z1_ll, z1_rr, -(gamma + alpha) / alpha)
  h5 = h1 * (
    gamma / (gamma - 1) * z5_power_mean * z1_power_expmean
    + 0.25 * (z2_ll + z2_rr)^2
    + 0.25 * (z3_ll + z3_rr)^2
    - 0.25 * (z2_ll^2 + z2_rr^2 + z3_ll^2 + z3_rr^2)
  )

  return SVector(h1, h2, h3, h5)
end
