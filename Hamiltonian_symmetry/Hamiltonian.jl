using LinearAlgebra, BenchmarkTools, Plots
gr()

const Nx, x₀, xₙ = 51, -50.0, 50.0
const Ny, y₀, yₙ = 51, -50.0, 50.0
const Nz, z₀, zₙ = 51, -50.0, 50.0

function RotMat(ϕ, u::Array{Float64, 1}=[0.0, 0.0, 1.0])
    ux, uy, uz = u
    [cos(ϕ)+ux^2*(1-cos(ϕ))     ux*uy*(1-cos(ϕ))-uz*sin(ϕ)  ux*uz*(1-cos(ϕ))+uy*sin(ϕ);
     uy*ux*(1-cos(ϕ))+uz*sin(ϕ) cos(ϕ)+uy^2*(1-cos(ϕ))      uy*uz*(1-cos(ϕ))-ux*sin(ϕ);
     uz*ux*(1-cos(ϕ))-uy*sin(ϕ) uz*uy*(1-cos(ϕ))+ux*sin(ϕ)  cos(ϕ)+uz^2*(1-cos(ϕ))
    ]
end

function RotMat(ϕ, α, β)
    RotMat(ϕ, axis(α, β))
end

function axis(ϕ, θ)
    [cos(ϕ)*sin(θ), sin(ϕ)*sin(θ), cos(θ)]
end

function Gaussian(x::Array{Float64, 1}, σ::Float64; A::Float64=-1.0, μ::Array{Float64, 1} = [0.,0.,0.])
    A/√(2π) * exp(-0.5 * (norm(x .- μ)/σ)^2)
end

function Piezo(x::Array{Float64, 1}, σ::Float64; A::Float64=-1.0, μ::Array{Float64, 1} = [0.,0.,0.], scale::Float64=0.5, rotmat)
    piezo = 0.0
    r = zeros(Float64, 3)
    for i in [-1.,1.], j in [-1.,1.], k in [-1.,1.]
        r[1] = i
        r[2] = j
        r[3] = k
        factor = -i*j*k
        piezo += Gaussian(x, σ*scale, A=A*factor, μ=μ .+ σ .* (rotmat * r))
    end
    return piezo
end

function Hamiltonian(σ::Float64, D::Float64; E_z::Float64=0.0, E_x::Float64=0.0, E_y::Float64=0.0, shift_x::Float64=0.0, shift_y::Float64=0.0, rotmat::Array{Float64, 2})
    E_vec = (rotmat * [E_x, E_y, E_z])'
    V = zeros(Float64, (Nx, Ny, Nz))
    r = zeros(Float64, 3)
    μ_bottom = rotmat * [-shift_x/2.0, -shift_y/2.0, -D/2.0]
    μ_top = rotmat * [shift_x/2.0, shift_y/2.0, D/2.0]
    for (i,x) in enumerate(range(x₀, xₙ, length=Nx)), (j,y) in enumerate(range(y₀, yₙ, length=Ny)), (k,z) in enumerate(range(z₀, zₙ, length=Nz))
        r .= [x,y,z]
        V[i, j, k] += Gaussian(r, σ, μ=μ_bottom) + Gaussian(r, σ, μ=μ_top)
        V[i, j, k] += Piezo(r, σ, μ=μ_bottom, rotmat=rotmat) + Piezo(r, σ, μ=μ_top, rotmat=rotmat)
        V[i, j, k] += -E_vec * r
    end
    return V
end

function measure(A::AbstractArray, B::AbstractArray)
    0.5 * my_norm(A .- B)^2 / (my_norm(A)^2 + my_norm(B)^2)
end

function my_norm(A::AbstractArray)
    √sum(A.^2)
end

# E_range = -1e-4:1e-5:1e-4
# result = [measure(Hamiltonian(20., 50., E_x = E, E_y = E, rotmat=RotMat(0, [0.,0.,1.])),
#                   Hamiltonian(20., 50., E_x = E, E_y = E, rotmat=RotMat(π, [0.,0.,1.]))) for E in E_range]
# plot(E_range, result)

ϕ_range = range(0.05, 0.2, length=21)
# ϕ_range = range(0, π/2, length=11)
E = 0e-4
Ez = 0e-3
shift = 10.0
result = [measure(Hamiltonian(20., 50., E_z = Ez, E_x = E, E_y = E, shift_x = shift, rotmat=RotMat(0, 0., ϕ)),
                  Hamiltonian(20., 50., E_z = Ez, E_x = E, E_y = E, shift_x = shift, rotmat=RotMat(π, 0., ϕ))) for ϕ in ϕ_range]
# result = [measure(Hamiltonian(20., 50., E_z = Ez, E_x = E, E_y = E, shift_x = shift/√2, shift_y = shift/√2, rotmat=RotMat(0, π/4, ϕ)),
#                   Hamiltonian(20., 50., E_z = Ez, E_x = E, E_y = E, shift_x = shift/√2, shift_y = shift/√2, rotmat=RotMat(π, π/4, ϕ))) for ϕ in ϕ_range]
print(minimum(result))
plot(ϕ_range, result, xlabel="θ (rad)", ylabel="Asymmetry")


# H₀ = Hamiltonian(20., 50., rotmat=RotMat(π/4, π/4, 0))
# heatmap(H₀[:,:,21])
# savefig("Changing axis.png")

# X = zeros(Nx, Nz)
# for i in 1:Nx
#     X[:, i] .= H₀[i, i, :]
# end
# heatmap(X)