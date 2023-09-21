# Definitions for Mie coeffcients of spheres (layered or not)
abstract type MultipoleOrder end

struct EM <: MultipoleOrder
    order :: Int
end
EM() = EM(1)


struct MM <: MultipoleOrder
    order :: Int
end
MM() = MM(1)


abstract type AbstractObject end

struct Sphere{T <: Real, U <: OptProp} <: AbstractObject
    material :: U
    radius   :: T
end

struct SphereCavity{T<:Real,U,V<:OptProp} <: AbstractObject
    material::U
    inner_material::V
    radius::T
end

permittivity(object::AbstractObject, w) = permittivity(object.material, w)
refractive_index(object::AbstractObject, w) = refractive_index(object.material, w)

# Helper functions
radius(object :: AbstractObject) = object.radius
size_parameter(object :: AbstractObject,w) = object.radius*w/c0

function max_order(object::AbstractObject, w)
    x = size_parameter(object, w)
    x_stop = x + 4 * x^(0.3333) + 2
    y_stop = abs(x * refractive_index(object, w))
    n_max = round(Int, max(x_stop, y_stop) + 15)
    return x_stop, y_stop, n_max
end

"""
    softinvsqrt(x::T) where {T}

Inverse square root.
"""
softinvsqrt(x::T) where {T} = sqrt(x) / (abs(x) + eps(float(T)))

"""
    softinvsqrt(x::Complex{T}) where {T}

Inverse square root for complex arguments.
"""
softinvsqrt(z::Complex{T}) where {T} = sqrt(z) / (abs(z) + eps(float(T)))

"""
    sphericalbesselj(nu, x)

Spherical bessel function of the first kind at order `nu`, ``j_ν(x)``. This is the non-singular
solution to the radial part of the Helmholz equation in spherical coordinates.
"""
sphericalbesselj(nu, x) = √(π/2) * softinvsqrt(x) * besselj(nu + 1//2, x)


"""
    sphericalbessely(nu, x)

Spherical bessel function of the second kind at order `nu`, ``y_ν(x)``. This is the singular
solution to the radial part of the Helmholz equation in spherical coordinates. Sometimes
known as a spherical Neumann function.
"""
sphericalbessely(nu, x) = √(π/2) * softinvsqrt(x) * bessely(nu + 1//2, x)


"""
    riccatibesselpsi(nu, x)

Riccati-Bessel function of the first kind at order `nu`, ``ψ_ν(x) = xj_ν(x)``.
"""
riccatibesselpsi(nu, x) = √(π*x/2) * besselj(nu + 1//2, x)

"""
    riccatibesselksi1(nu, x)

Riccati-Bessel function of the third kind at order `nu`, ``ξ_ν(x) = x(j_ν(x) + iy_ν(x))``.
"""
riccatibesselksi1(nu, x) = √(π*x/2) * (besselj(nu + 1//2, x) + im * bessely(nu + 1//2, x))

function mie_coefficients_openlibm(sphere::Sphere, w)
    x_stop, y_stop, n_max = max_order(sphere, w)
    size_par = size_parameter(sphere, w)
    ref_idx = refractive_index(sphere, w)
    y = size_par * ref_idx
    T = typeof(size_par)

    d = zeros(Complex{T}, n_max + 1)
    av = zeros(Complex{T}, round(Int, x_stop))
    bv = zeros(Complex{T}, round(Int, x_stop))

    # downward recursion
    @inbounds for n = n_max-1:-1:1
        rn = n + 1
        d[n] = (rn / y) - (1 / (d[n+1] + rn / y))
    end

    # upward recursion
    @inbounds for n = 1:round(Int, x_stop)
        psi_n = riccatibesselpsi(n, size_par)
        ksi_n = riccatibesselksi1(n, size_par)
        psi_n_1 = riccatibesselpsi(n - 1, size_par)
        ksi_n_1 = riccatibesselksi1(n - 1, size_par)
        t_a = d[n] / ref_idx + n / size_par
        t_b = d[n] * ref_idx + n / size_par
        av[n] = (t_a * psi_n - psi_n_1) / (t_a * ksi_n - ksi_n_1)
        bv[n] = (t_b * psi_n - psi_n_1) / (t_b * ksi_n - ksi_n_1)
    end
    return av, bv
end


function mie_coefficients_handcoded(sphere::Sphere, w)
    x_stop, y_stop, n_max = max_order(sphere, w)
    size_par = size_parameter(sphere, w)
    ref_idx = refractive_index(sphere, w)
    y = size_par * ref_idx

    d = zeros(Complex{Float64}, n_max + 1)
    av = zeros(Complex{Float64}, round(Int, x_stop))
    bv = zeros(Complex{Float64}, round(Int, x_stop))

    # downward recursion
    @inbounds for n = n_max-1:-1:1
        rn = n + 1
        d[n] = (rn / y) - (1 / (d[n+1] + rn / y))
    end

    psi0 = cos(size_par)
    psi1 = sin(size_par)
    chi0 = -sin(size_par)
    chi1 = cos(size_par)

    xi1 = Complex{Float64}(psi1, -chi1)

   # upward recursion
    @inbounds for n = 1:round(Int, x_stop)
        psi = (2 * n - 1) * psi1 / size_par - psi0
        chi = (2 * n - 1) * chi1 / size_par - chi0
        xi = Complex{Float64}(psi, -chi)
        t_a = d[n] / ref_idx + n / size_par
        t_b = d[n] * ref_idx + n / size_par
        av[n] = (t_a * psi - psi1) / (t_a * xi - xi1)
        bv[n] = (t_b * psi - psi1) / (t_b * xi - xi1)

        psi0 = psi1
        psi1 = psi
        chi0 = chi1
        chi1 = chi
        xi1 = Complex{Float64}(psi1, -chi1)
    end

    return av, bv
end

"""
    mie_coefficients(sphere::Sphere, w; handcoded=false)

Compute mie coefficients ``a_n`` and ``b_n`` for a sphere at frequency `w`.

If `handcoded` is `false`, it uses the bessel functions from SpecialFunctions.
Otherwise it uses the code from Bohren and Huffman which is about 10 times
but not reliable for orders `nu = around`.
"""
function mie_coefficients(sphere::Sphere, w; handcoded=false )
    if handcoded==false
        return  mie_coefficients_openlibm(sphere::Sphere, w)
    elseif handcoded==true
        return  mie_coefficients_handcoded(sphere::Sphere, w)
    end
end



"""
    efficiencies(sphere::Sphere, w; hardcoded=true)

Compute scattering, extinction and absorption efficiencies
for a sphere at frequency `w`.

If `hardcoded` is `true`, Mie coefficients are obtained using the code from
Bohren and Huffman, otherwise it uses a call to the bessel functions from SpecialFunctions.
The hardcoded version is about 10 times faster.
"""
function efficiencies(sphere::Sphere, w; kwargs...)
    size_par = size_parameter(sphere, w)
    x_stop, y_stop, n_max = max_order(sphere, w)

    a, b = mie_coefficients(sphere, w; kwargs...)
    q_sca = 0
    q_ext = 0

    for n = 1:round(Int, x_stop)
        q_sca += (2 * n + 1) * (abs(a[n])^2 + abs(b[n])^2)
        q_ext += (2 * n + 1) * real(a[n] + b[n])
    end

    return 2 / (size_par^2) * q_sca,
    2 / (size_par^2) * q_ext,
    2 / (size_par^2) * (q_ext - q_sca)
end


function efficiencies(sphere::Sphere, comp :: MultipoleOrder  , w; kwargs...)
    size_par = size_parameter(sphere, w)
    x_stop, y_stop, n_max = max_order(sphere, w)

    a, b = mie_coefficients(sphere, w; kwargs...)
    q_sca = 0
    q_ext = 0

    if typeof(comp) == EM
        c = a
    elseif typeof(comp) == MM
        c = b
    end

    for n = 1:round(Int, x_stop)
        q_sca += (2 * n + 1) * (abs(c[n])^2)
        q_ext += (2 * n + 1) * real(c[n])
    end

    return 2 / (size_par^2) * q_sca,
    2 / (size_par^2) * q_ext,
    2 / (size_par^2) * (q_ext - q_sca)
end
