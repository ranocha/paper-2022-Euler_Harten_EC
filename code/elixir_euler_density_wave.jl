
using LinearAlgebra
using Plots
using Trixi

function compute_spectrum(volume_flux, alpha)
  # semidiscretization of the compressible Euler equations
  equations = HartensCompressibleEulerEquations2D(1.4, alpha)

  initial_condition = initial_condition_density_wave

  solver = DGSEM(polydeg=5, surface_flux=volume_flux,
                 volume_integral=VolumeIntegralFluxDifferencing(volume_flux))

  coordinates_min = (-1.0, -1.0)
  coordinates_max = ( 1.0,  1.0)
  mesh = TreeMesh(coordinates_min, coordinates_max,
                  initial_refinement_level=2,
                  n_cells_max=30_000)

  semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver)

  # compute the spectrum
  J = jacobian_ad_forward(semi)
  λ = eigvals(J)
end

function plot_spectrum(volume_flux, alpha, filename)
  λ = compute_spectrum(volume_flux, alpha)

  fontsizes = (xtickfontsize=18, ytickfontsize=18,
               xguidefontsize=20,  yguidefontsize=20, legendfontsize=18)
  kwargs = (xguide="\$\\mathrm{Re}(\\lambda)\$", yguide="\$\\mathrm{Im}(\\lambda)\$",
            legend=:none, size=(700, 500), fontsizes...)
  scatter(real.(λ), imag.(λ); kwargs...)
  savefig(filename)
end

# β = (α + γ) / (1 - γ) = 1 for α = 1 - 2γ (γ = 1.4)
plot_spectrum(flux_sjögreen_yee, -1.8,
              joinpath(dirname(@__DIR__), "figures", "spectrum_beta_1.pdf"))

# β = (α + γ) / (1 - γ) = 2 for α = 2 - 3γ (γ = 1.4)
plot_spectrum(flux_sjögreen_yee, -2.2,
              joinpath(dirname(@__DIR__), "figures", "spectrum_beta_2.pdf"))

# β = (α + γ) / (1 - γ) = -6 for α = 1 (γ = 1.4)
plot_spectrum(flux_sjögreen_yee, 1.0,
              joinpath(dirname(@__DIR__), "figures", "spectrum_beta_m6.pdf"))
