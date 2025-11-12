## Density estimation example
# from "Non-Parametric Models for Non-Negative Functions"
#
# Lorenzo Shaikewitz

using LinearAlgebra
using Random
using Plots

using Convex
using MosekTools

# parameters
λ = 1e-3 # regularization
σ = 1. # kernel variance
s = 1. # kernel scale
n = 50 # number of samples

# distribution parameters
# use a mixture of gaussians (μ1, σ1) and (μ2, σ2) mixed with prob. p.
μ1 = -1.
μ2 = 1.
σ2_1 = 0.3^2
σ2_2 = 0.3^2
p = 0.5

# base kernel (Gaussian)
k(x, y) = s*exp(-((x-y)'*(x-y)) / (2σ^2))

Random.seed!(0)
# generate independent samples from Gaussian mixture
function rand_gmm(n)
    indicator = rand(n)
    out  = (indicator .> p).*(sqrt(σ2_1)*randn(n) .+ μ1)
    out += (indicator .≤ p).*(sqrt(σ2_2)*randn(n) .+ μ2)
end
# this is our "calibration" data (iid samples from distribution)
x_cal = rand_gmm(n)

## Optimization problem
# the goal is to estimate P(x | x_cal)
# we'll use a semidefinite model
B = Semidefinite(n)
v(x) = k.([x], x_cal)
# Model: P(x | x_cal) = v(x)'*B*v(x)

# maximum likelihood estimate
# max P(x | x_cal)
# take the -log likelihood:
objective  = -1/n * sum([log(v(xi_cal)'*B*v(xi_cal)) for xi_cal in x_cal])
# regularization term
objective += λ*(nuclearnorm(B) + 0.005*tr(B))

# integral constraint
# enforces integral of P(x | x_cal) over x = 1
# M is just the integral of the product of kernels
M = s^2*exp.(-1/(4σ^2)*((x_cal .- x_cal').^2))*sqrt(σ^2*π) # closed form from kernel
constraints = [tr(B*M) == 1]

# solve!
problem = minimize(objective, constraints)
solve!(problem, Mosek.Optimizer)

# we have our model!
B_val = evaluate(B)

# ## Mode estimation
# x = Variable()
# obj_mode = sum([-(1/(2σ^2))*((x-x_cal[i])^2 + (x-x_cal[j])^2) for i = 1:n, j = 1:n])
# prob_mode = maximize(obj_mode)
# solve!(prob_mode, Mosek.Optimizer)


## plot
pdf_n(x, μ, σ2) = 1/(sqrt(2π*σ2))*exp(-(x-μ)^2 / (2σ2))
pdf(x) = p*pdf_n(x, μ1, σ2_1) + (1-p)*pdf_n(x, μ2, σ2_2)

x_dense = collect(-4:0.01:4)
Plots.plot(x_dense, transpose.(v.(x_dense)) .* [B_val] .* v.(x_dense), label="learned",lw=2)
Plots.plot!(x_dense, pdf.(x_dense), label="ground truth", ylabel="p(x | cal)",lw=2)
Plots.scatter!(x_cal, 0 .*x_cal, label="training points", markershape=:x, msw=2)