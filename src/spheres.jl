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
    temperature :: T
end

struct Cavity{T <: Real, U <: OptProp, V <: OptProp} <: AbstractObject
    material :: U
    inner_material :: V
    radius   :: T
    temperature :: T
end


permittivity(object::AbstractObject, w) = permittivity(object.material, w)
refractive_index(object::AbstractObject, w) = refractive_index(object.material, w)

# Helper functions
radius(object :: AbstractObject) = object.radius
size_parameter(object :: AbstractObject,w) = object.radius*w/c0
temperature(object::AbstractObject) = object.temperature

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


"""
    riccatibesselksi2(nu, x)

Riccati-Bessel function of the third kind at order `nu`, ``ξ_ν(x) = x(j_ν(x) - iy_ν(x))``.
"""
riccatibesselksi2(nu, x) = √(π*x/2) * (besselj(nu + 1//2, x) - im * bessely(nu + 1//2, x))


function mie_coefficients_openlibm(sphere::Sphere, w, nmax = 0)
    x_stop, y_stop, n_max = max_order(sphere, w)
    nmax != 0 ? n_max = nmax : nothing
    size_par = size_parameter(sphere, w)
    ref_idx = refractive_index(sphere, w)
    y = size_par * ref_idx
    T = typeof(size_par)

    d = zeros(Complex{T}, n_max + 1)
    av = zeros(Complex{T}, n_max)
    bv = zeros(Complex{T}, n_max)

    # downward recursion
    @inbounds for n = n_max:-1:1
        rn = n + 1
        d[n] = (rn / y) - (1 / (d[n+1] + rn / y))
    end

    # upward recursion
    @inbounds for n = 1:n_max #round(Int, x_stop)
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


function mie_coefficients_handcoded(sphere::Sphere, w, nmax = 0 )
    x_stop, y_stop, n_max = max_order(sphere, w)
    nmax != 0 ? n_max = nmax : nothing
    size_par = size_parameter(sphere, w)
    ref_idx = refractive_index(sphere, w)
    y = size_par * ref_idx

    d = zeros(Complex{Float64}, n_max + 1)
    av = zeros(Complex{Float64}, round(Int, x_stop))
    bv = zeros(Complex{Float64}, round(Int, x_stop))

    # downward recursion
    @inbounds for n = n_max-1:-1:1
        rn = n + 1
        d[n] = (rn / y) - (1 / (d[rn] + rn / y))
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
function mie_coefficients(sphere::Sphere, w, nmax = 0; handcoded=false )
    if handcoded==false
        return  mie_coefficients_openlibm(sphere::Sphere, w, nmax)
    elseif handcoded==true
        return  mie_coefficients_handcoded(sphere::Sphere, w, nmax)
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




function mie_coefficients(cavity::Cavity, w, nmax = 0)
    x_stop, y_stop, n_max = max_order(cavity, w)
    nmax != 0 ? n_max = nmax : nothing
    x = size_parameter(cavity, w)
    ref_idx = refractive_index(cavity, w)
    y = x * ref_idx
    T = typeof(x)

    gx = zeros(Complex{T}, n_max + 1)
    gx[1] = -im
    gy = zeros(Complex{T}, n_max + 1)
    gy[1] = -im
    dx = zeros(Complex{T}, n_max + 1)
    cv = zeros(Complex{T}, n_max)
    dv = zeros(Complex{T}, n_max)
    

    # downward recursion
    @inbounds for n = n_max:-1:1
        rn = n + 1
        dx[n] = (rn / x) - (1 / (dx[rn] + rn / x))
    end
        
    # upward recursion
    @inbounds for n = 1:n_max
        if n==1
            gy[1] = im # check this sign!!! Wy originally -im??? 
        else    
            rn = n - 1 
            gy[n] = (1 / (- gy[rn] + n / y)) - (n / y)
        end

        psi_n = riccatibesselpsi(n, x)
        ksi_n = riccatibesselksi1(n, x)
        psi_n_1 = riccatibesselpsi(n - 1, x)
        ksi_n_1 = riccatibesselksi1(n - 1, x)
          
        t_c = gy[n] + ref_idx * n / x
        t_d = gy[n] * ref_idx + n / x

        
        cv[n] = -( t_c * ksi_n - ref_idx * ksi_n_1 ) / ( t_c * psi_n - ref_idx * psi_n_1 ) #TE
        dv[n] = -( t_d * ksi_n - ksi_n_1 ) / ( t_d * psi_n - psi_n_1 ) # TM
    end
    return cv, dv
end


struct SphericalShell{U,V} <: AbstractObject
    sphere::U
    cavity::V
end

sphere(spsh :: SphericalShell) = spsh.sphere
cavity(spsh :: SphericalShell) = spsh.cavity

function transmission(sphericalshell :: SphericalShell, w; nmax=0, rtol=1e-6)
    (;sphere, cavity) = sphericalshell
    xsph = size_parameter(sphere, w)
    xsph_stop, ysph_stop, nsph_max = max_order(sphere, w)
    
    xcav = size_parameter(cavity, w)
    xcav_stop, ycav_stop, ncav_max = max_order(cavity, w)
   
    n_max = min(nsph_max, ncav_max)
    nmax != 0 ? n_max = nmax : nothing

    te_sph, tm_sph = mie_coefficients(sphere, w, n_max)
    te_cav, tm_cav = mie_coefficients(cavity, w, n_max)
    
    qw_te = 0.0
    qw_tm = 0.0
    
    qw_te_1 = 0.0
    qw_tm_1 = 0.0
    err = 1.0

    if nmax != 0
        for i in 1:nmax
            #qw_te += (2*i+1) * (real(te_cav[i]) + 1) * (real(te_sph[i]) - abs(te_sph[i])^2) / abs(1 - te_cav[i]*te_sph[i])^2
            #qw_tm += (2*i+1) * (real(tm_cav[i]) + 1) * (real(tm_sph[i]) - abs(tm_sph[i])^2) / abs(1 - tm_cav[i]*tm_sph[i])^2       
            qw_te += (2*i+1) * exp(log(real(te_cav[i]) + 1) + log(real(te_sph[i]) - abs(te_sph[i])^2) - 2*log(abs(1 - exp(log(te_cav[i]) + log(te_sph[i]))))) #check signs!!
            qw_tm += (2*i+1) * exp(log(real(tm_cav[i]) + 1) + log(real(tm_sph[i]) - abs(tm_sph[i])^2) - 2*log(abs(1 - exp(log(tm_cav[i]) + log(tm_sph[i])))))        
        end
        return (te =  qw_te, tm = qw_tm, total = (qw_te + qw_tm), n_max = n_max)
    else
        i=0
        while (err > rtol) && (i < n_max)
            i += 1 
            qw_te += (2*i+1) * exp(log(real(te_cav[i]) + 1) + log(real(te_sph[i]) - abs(te_sph[i])^2) - 2*log(abs(1 - exp(log(te_cav[i]) + log(te_sph[i]))))) #check signs!!
            qw_tm += (2*i+1) * exp(log(real(tm_cav[i]) + 1) + log(real(tm_sph[i]) - abs(tm_sph[i])^2) - 2*log(abs(1 - exp(log(tm_cav[i]) + log(tm_sph[i])))))        
            err = abs(qw_te_1 + qw_tm_1 - qw_te - qw_tm)/abs(qw_te + qw_tm)
            qw_te_1 = qw_te
            qw_tm_1 = qw_tm 
        end
        i > n_max ? error("The calculation did not converge. The number of iterations exceeds the value of n_max.") : nothing 
        err > rtol ? error("The calculation did not converge. The relative error is larger than rtol.") : nothing
        return (te =  qw_te, tm = qw_tm, total = (qw_te + qw_tm), numiter = i, n_max = n_max, rel_error = err)
    end
end


"Spectral heat flux"
function heat_flux(sphericalshell :: SphericalShell, w; nmax=0, rtol=1e-6)
    (;sphere, cavity) = sphericalshell
    T1 = temperature(sphere)
    T2 = temperature(cavity)
    delTheta = (bose_einstein(w,T1) - bose_einstein(w,T2))
    transm = transmission(sphericalshell, w; nmax, rtol)
    qw_te = delTheta * transm[:te]
    qw_tm = delTheta * transm[:tm]
    qw_tot = delTheta * transm[:total]
    return (te = qw_te, tm = qw_tm, total = qw_tot, nmax = transm[:n_max])
end


"Total heat flux"
function heat_flux(sphericalshell :: SphericalShell; w1 = 0.0, w2 = Inf, nmax=0, rtol=1e-6, rtol_int = 1e-6)
    (;sphere, cavity) = sphericalshell
    hf_w(w) = heat_flux(sphericalshell, w; nmax, rtol)
    qte_w(u) = hf_w(u*kb/ħ)[:te]
    qtm_w(u) = hf_w(u*kb/ħ)[:tm]
    n_max_w = hf_w(w2)[:nmax]
    
    u1 = w1*ħ/kb
    u2 = w2*ħ/kb

    (qte, err) = quadgk(qte_w, u1 , u2 ; rtol = rtol_int)
    (qtm, err) = quadgk(qtm_w, u1 , u2 ; rtol = rtol_int)
    return (te = qte * kb / ħ, tm = qtm * kb / ħ, total = (qte + qtm) * kb / ħ, nmax = n_max_w)
end