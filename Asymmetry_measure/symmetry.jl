using LinearAlgebra
using SparseArrays
using PETScBinaryIO

const h = 40
const hp1 = 41

# const hz = 40
# const hzp1 = 41
const hz = 50
const hzp1 = 51

# f0001 - bad crop
# f0002 - no crop
# f0004 - good crop
# const filename = "results/elong/Hamiltonians/s0001f0004.dat" # Symmetrical
# const filename = "results/elong/Hamiltonians/s0002f0004.dat" # Elong110
# const filename = "results/elong/Hamiltonians/s0003f0004.dat" # Elong100

function i_to_coord(i)
    i = i-1 # Start at 0
    kp = i % 8 + 1
    i = ceil(Int, (i+1)/8) - 1
    x = i % h + 1
    i = ceil(Int, (i+1)/h) - 1
    y = i % h + 1
    z = ceil(Int, (i+1)/h)

    return (kp, x, y, z)
end

function coord_to_i(kp, x, y, z)
    return (((z-1)*h+(y-1))*h+(x-1))*8+kp
end

function flip_110(kp, x, y, z)
    return (kp, y, x, z)
end

function flip_m110(kp, x, y, z)
    return (kp, hp1-y, hp1-x, z)
end

function flip_100(kp, x, y, z)
    return (kp, x, hp1-y, z)
end

function flip_010(kp, x, y, z)
    return (kp, hp1-x, y, z)
end

function rotate_180(kp, x, y, z)
    return (kp, hp1-x, hp1-y, z)
end

function symmetry_metric(A,B)
    return 0.5 * norm(A - B) / (norm(A) + norm(B))
end

function Jx()
    result = zeros(ComplexF64, 8, 8)
    result[1,2] = 1
    result[2,1] = 1

    result[3,4] = 3^0.5
    result[4,3] = 3^0.5
    result[4,5] = 2
    result[5,4] = 2
    result[5,6] = 3^0.5
    result[6,5] = 3^0.5

    result[7,8] = 1
    result[8,7] = 1

    return 0.5 * result
end

function Jy()
    result = zeros(ComplexF64, 8, 8)
    result[1,2] = -1
    result[2,1] = 1

    result[3,4] = -3^0.5
    result[4,3] = 3^0.5
    result[4,5] = -2
    result[5,4] = 2
    result[5,6] = -3^0.5
    result[6,5] = 3^0.5

    result[7,8] = -1
    result[8,7] = 1

    return 0.5im * result
end

function Jz()
    result = zeros(ComplexF64, 8, 8)
    result[1,1] = 0.5
    result[2,2] = -0.5
    result[3,3] = 1.5
    result[4,4] = 0.5
    result[5,5] = -0.5
    result[6,6] = -1.5
    result[7,7] = 0.5
    result[8,8] = -0.5
    return result
end

function inversion_orb()
    result = spzeros(ComplexF64, 8, 8)
    result[1, 1] = 1.0
    result[2, 2] = 1.0

    result[3, 3] = -1.0
    result[4, 4] = -1.0
    result[5, 5] = -1.0
    result[6, 6] = -1.0

    result[7, 7] = -1.0
    result[8, 8] = -1.0
    return result
end

function rotate(axis)
    if isa(axis, Number)
        axis = deg2rad(axis)
        J = Jx() * cos(axis) + Jy() * sin(axis)
    elseif axis == "110"
        J = (Jx() + Jy()) * 2^(-0.5)
    elseif  axis == "m110"
        J = (-Jx() + Jy()) * 2^(-0.5)
    elseif  axis == "100"
        J = Jx()
    elseif  axis == "010"
        J = Jy()
    elseif  axis == "001"
        J = Jz()
    else
        J = spzeros(ComplexF64, 8, 8) # Will crash if applied
    end
    X = sparse(exp(-1im*pi*J))
    for (x, y, v) in zip(findnz(X)...)
        if abs(v) < 1e-10
            X[x, y] = 0
        end
    end
    dropzeros!(X)
    return X
end

function kp_8x8(axis)
    return rotate(axis) * inversion_orb()
end

function generate_A()    
    A_110 = spzeros(ComplexF64, h^2*hz, h^2*hz)
    A_m110 = spzeros(ComplexF64, h^2*hz, h^2*hz)
    A_rotate = spzeros(ComplexF64, h^2*hz, h^2*hz)

    for i in 1:8:8*h^2*hz
        A_110[coord_to_i(flip_110(i_to_coord(i)...)...)÷8+1, i÷8+1] = 1.0
        A_m110[coord_to_i(flip_m110(i_to_coord(i)...)...)÷8+1, i÷8+1] = 1.0
        A_rotate[coord_to_i(rotate_180(i_to_coord(i)...)...)÷8+1, i÷8+1] = 1.0
    end

    A_110 = kron(A_110, kp_8x8("m110"))
    A_m110 = kron(A_m110, kp_8x8("110"))
    A_rotate = kron(A_rotate, rotate("001"))
    
    println("A generation done")
    return A_110, A_m110, A_rotate
end

function generate_A(alpha)    
    A_100 = spzeros(ComplexF64, h^2*hz, h^2*hz)
    A_010 = spzeros(ComplexF64, h^2*hz, h^2*hz)
    A_rotate = spzeros(ComplexF64, h^2*hz, h^2*hz)

    for i in 1:8:8*h^2*hz
        A_100[coord_to_i(flip_100(i_to_coord(i)...)...)÷8+1, i÷8+1] = 1.0
        A_010[coord_to_i(flip_010(i_to_coord(i)...)...)÷8+1, i÷8+1] = 1.0
        A_rotate[coord_to_i(rotate_180(i_to_coord(i)...)...)÷8+1, i÷8+1] = 1.0
    end

    A_100 = kron(A_100, kp_8x8(90 + alpha))
    A_010 = kron(A_010, kp_8x8(alpha))
    A_rotate = kron(A_rotate, rotate("001"))
    
    println("A generation done")
    return A_100, A_010, A_rotate
end

function show_row(H, i)
    for (x,v) in zip(findnz(H[i,:])...)
        y = i_to_coord(x)
        println(x, y, v)
    end    
end

function print_diagonal(H)
    for (x,v) in zip(findnz(diag(H))...)
        y = i_to_coord(x)
        println(x, y, v)
    end
end

function test_coord_change()
    x = coord_to_i(1, 1, 1, 2)
    x = coord_to_i(1, 40, 40, 40)
    @show x
    y = i_to_coord(x)
    @show y

    x = (1, 1, 1, 1)
    @show x
    @show rotate_180(x...)
    @show flip_110(x...)
    @show flip_m110(x...)    
end

function test_matrices()
    @show dropzeros(sparse(Jx()))
    @show dropzeros(sparse(Jy()))
    @show rotate("110")
    @show rotate("m110")
    @show kp_8x8("110")
    @show kp_8x8("m110")
    @show kp_8x8("110") * kp_8x8("m110")
end

function read_H(filename)
    H = readpetsc(filename, int_type=Int32, scalar_type=ComplexF64)
    H = H[1]
    dropzeros!(H)
    println("Hamiltonian read done")
    return H
end

function show_elements(H, H_110, H_m110, H_rotate, kp, x, y, z)
    show_row(H, coord_to_i(kp, x, y, z))
    show_row(H_110, coord_to_i(kp, x, y, z))
    show_row(H_m110, coord_to_i(kp, x, y, z))
    show_row(H_rotate, coord_to_i(kp, x, y, z))
end

function calc_metrics(filename, A_110, A_m110, A_rotate)
    @show filename
    H = read_H(filename)

    # H = abs.(H)

    H_110 = A_110 * H * A_110'
    H_m110 = A_m110 * H * A_m110'
    H_rotate = A_rotate * H * A_rotate'

    # show_elements(H, H_110, H_m110, H_rotate, 1, 40, 40, 50)

    metric_110 = symmetry_metric(H, H_110)
    metric_m110 = symmetry_metric(H, H_m110)
    metric_rotate =  symmetry_metric(H, H_rotate)
    @show metric_110
    @show metric_m110
    @show metric_rotate    

    return (metric_110, metric_m110, metric_rotate)
end

# A_110, A_m110, A_rotate = generate_A()

# calc_metrics("results/elong/Hamiltonians/s0002f0004.dat", A_110, A_m110, A_rotate)
# calc_metrics("results/rotate_low/Hamiltonians/s0013f0001.dat", A_110, A_m110, A_rotate)
# calc_metrics("results/rotate_low/Hamiltonians/s0013f0001.dat", A_100, A_010, A_rotate)

for i in 1:14
    # A_100, A_010, A_rotate = generate_A()
    A_100, A_010, A_rotate = generate_A(i<=10 ? 5*(i-1) : 30+i)
    # filename = "results/elong_100_low/Hamiltonians/s00"
    # filename = "results/rotate_low/Hamiltonians/s00"
    filename = "results/rotate_shift/Hamiltonians/s00"
    if i<10
        filename = filename * "0"
    end
    filename = filename * string(i) * "f0001.dat"
    # filename = filename * string(i) * "f0002.dat"
    x, y, z = calc_metrics(filename, A_100, A_010, A_rotate)
    open("Hsym_results_cut.txt", "a") do io
        write(io, string(i) * "\t" * string(x) * "\t" * string(y) * "\t" * string(z) * "\n")
    end
end