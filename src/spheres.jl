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


a(n, nu, x)= (-1)^(n+1)*2*(nu + n - 1.0)/x


function cont_frac(n, nu, x, stop=1)
    if n>stop+1
       return a(n, nu, x) + 1/cont_frac(n-1, nu, x, stop)
    else n==stop+1
       return a(n, nu, x) + 1/a(n-1, nu, x)
    end
 end


 function start_value(n, x)
    i = 2
    a1 = a(1, n + 0.5, x)
    a2 = a(2, n + 0.5, x)
    a21 = a2 + 1.0/a1
    bessel_ratio = a1*a21/a2
    while a1 != a2 
        i += 1
        a1 = cont_frac(i, n + 0.5, x, 1)
        a2 = cont_frac(i, n + 0.5, x, 2)
        bessel_ratio *= a1/a2
    end
    return  -n/x + bessel_ratio
 end



function mie_coefficients_openlibm(sphere::Sphere, w, nmax = 0)
    x_stop, y_stop, n_max = max_order(sphere, w)
    n_max_down = n_max + 1
    nmax != 0 ? n_max = nmax : nothing
    nmax != 0 ? n_max_down = nmax : nothing
    size_par = size_parameter(sphere, w)
    ref_idx = refractive_index(sphere, w)
    y = size_par * ref_idx
    T = typeof(size_par)

    d = zeros(Complex{T}, n_max_down)
    d[n_max_down] = start_value(n_max_down, y)
    av = zeros(Complex{T}, n_max)
    bv = zeros(Complex{T}, n_max)
    
    # downward recursion
    for n = (n_max_down-1):-1:1
        rn = n + 1
        d[n] = (rn / y) - (1 / (d[n+1] + rn / y))
    end

    # upward recursion
    for n = 1:n_max #round(Int, x_stop)
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



function log_mie_coefficients(sphere::Sphere, w, nmax = 0)
    x_stop, y_stop, n_max = max_order(sphere, w)
    n_max_down = n_max + 1
    nmax != 0 ? n_max = nmax : nothing
    nmax != 0 ? n_max_down = nmax : nothing
    size_par = size_parameter(sphere, w)
    ref_idx = refractive_index(sphere, w)
    y = size_par * ref_idx
    T = typeof(size_par)

    d = zeros(Complex{T}, n_max_down)
    d[n_max_down] = start_value(n_max_down, y)
    logav = zeros(Complex{T}, n_max)
    logbv = zeros(Complex{T}, n_max)
    logconjav = zeros(Complex{T}, n_max)
    logconjbv = zeros(Complex{T}, n_max)
    
    # downward recursion
    for n = (n_max_down-1):-1:1
        rn = n + 1
        d[n] = (rn / y) - (1 / (d[n+1] + rn / y))
    end

    # upward recursion
    for n = 1:n_max #round(Int, x_stop)
        psi_n = riccatibesselpsi(n, size_par)
        ksi_n = riccatibesselksi1(n, size_par)
        psi_n_1 = riccatibesselpsi(n - 1, size_par)
        ksi_n_1 = riccatibesselksi1(n - 1, size_par)
        t_a = d[n] / ref_idx + n / size_par
        t_b = d[n] * ref_idx + n / size_par
        logav[n] = log(t_a * psi_n - psi_n_1) - log(t_a * ksi_n - ksi_n_1)
        logbv[n] = log(t_b * psi_n - psi_n_1) - log(t_b * ksi_n - ksi_n_1)
        logconjav[n] = log(conj(t_a * psi_n - psi_n_1)) - log(conj(t_a * ksi_n - ksi_n_1))
        logconjbv[n] = log(conj(t_b * psi_n - psi_n_1)) - log(conj(t_b * ksi_n - ksi_n_1))

    end
    return logav, logbv, logconjav, logconjbv
end




function mie_coefficients_handcoded(sphere::Sphere, w, nmax = 0 )
    x_stop, y_stop, n_max = max_order(sphere, w)
    n_max_down = n_max + 1
    nmax != 0 ? n_max = nmax : nothing
    nmax != 0 ? n_max_down = nmax : nothing
    size_par = size_parameter(sphere, w)
    ref_idx = refractive_index(sphere, w)
    y = size_par * ref_idx

    d = zeros(Complex{Float64}, n_max_down)
    d[n_max_down] = start_value(n_max_down, y)
    av = zeros(Complex{Float64}, n_max)
    bv = zeros(Complex{Float64}, n_max)

    # downward recursion
    for n = (n_max_down-1):-1:1
        rn = n + 1
        d[n] = (rn / y) - (1 / (d[rn] + rn / y))
    end

    psi0 = cos(size_par)
    psi1 = sin(size_par)
    chi0 = -sin(size_par)
    chi1 = cos(size_par)

    xi1 = Complex{Float64}(psi1, -chi1)

   # upward recursion
    for n = 1:n_max
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
    
    isnan(q_ext) ? (@show w) : nothing
    isnan(q_sca) ? (@show w) : nothing
    
    return 2 / (size_par^2) * q_sca,
    2 / (size_par^2) * q_ext,
    2 / (size_par^2) * (q_ext - q_sca)
end


function efficiencies_n(sphere::Sphere, w, order; kwargs...)
    size_par = size_parameter(sphere, w)
    
    a, b = mie_coefficients(sphere, w, order; kwargs...)
    q_sca = (2 * order + 1) * (abs(a[order])^2 + abs(b[order])^2)
    q_ext = (2 * order + 1) * real(a[order] + b[order])
    
    return 2 / (size_par^2) * q_sca,
    2 / (size_par^2) * q_ext,
    2 / (size_par^2) * (q_ext - q_sca)
end


function efficiencies_n(sphere::Sphere, w, mp_order :: MultipoleOrder; kwargs...)
    size_par = size_parameter(sphere, w)
    x_stop, y_stop, n_max = max_order(sphere, w)

    a, b = mie_coefficients(sphere, w; kwargs...)
    if typeof(mp_order) == EM
        c = a
    elseif typeof(mp_order) == MM
        c = b
    end
    q_sca = (2 * n + 1) * (abs(c[n])^2)
    q_ext = (2 * n + 1) * real(c[n])

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

    gy = zeros(Complex{T}, n_max + 1)
    dx = zeros(Complex{T}, n_max + 1)
    cv = zeros(Complex{T}, n_max)
    dv = zeros(Complex{T}, n_max)


    # downward recursion
    @inbounds for n = n_max:-1:1
        rn = n + 1
        dx[n] = (rn / x) - (1 / (dx[rn] + rn / x))
    end
    
    gy0 = im #riccatibesselksi1(-1, y)/riccatibesselksi1(0, y)

    # upward recursion
    @inbounds for n = 1:nmax
        if n == 1 
            gy[n] = 1/(- gy0 + n / y ) - n/y
        else
            rn = n - 1 
            gy[n] = 1/(- gy[rn] + n / y ) - n/y
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
    return  cv, dv
end


function analytical_mie_coefficients(cavity::Cavity, w, nmax)
    x = size_parameter(cavity, w)
    ref_idx = refractive_index(cavity, w)
    y = x * ref_idx
    T = typeof(x)

    cv_an = zeros(Complex{T}, n_max)
    dv_an = zeros(Complex{T}, n_max)
    
  
    # upward recursion
    for n in 1:nmax
        if n == 1 
            gy[n] = 1/(- gy0 + n / y ) - n/y
        else
            rn = n - 1 
            gy[n] = 1/(- gy[rn] + n / y ) - n/y
        end

        psi_n = riccatibesselpsi(n, x)
        ksi_n = riccatibesselksi1(n, x)
        ksiy_n = riccatibesselksi1(n, y)
        psi_n_1 = riccatibesselpsi(n - 1, x)
        ksi_n_1 = riccatibesselksi1(n - 1, x)
        ksiy_n_1 = riccatibesselksi1(n - 1, y)
        
        psip_n = psi_n_1 - n/x * psi_n
        ksip_n = ksi_n_1 - n/x * ksi_n
        ksipy_n = ksiy_n_1 - n/y * ksiy_n
       
        cv_an[n] = - (ksi_n * ksipy_n - ref_idx * ksip_n * ksiy_n) / (psi_n * ksipy_n - ref_idx * psip_n * ksiy_n) #TE
        dv_an[n] = - (ref_idx * ksi_n * ksipy_n - ksip_n * ksiy_n) / (ref_idx * psi_n * ksipy_n - psip_n * ksiy_n) #TM
    end
    return  cv_an, dv_an
end


function log_mie_coefficients(cavity::Cavity, w, nmax = 0)
    x_stop, y_stop, n_max = max_order(cavity, w)
    nmax != 0 ? n_max = nmax : nothing
    x = size_parameter(cavity, w)
    ref_idx = refractive_index(cavity, w)
    y = x * ref_idx
    T = typeof(x)

    gy = zeros(Complex{T}, n_max + 1)
    dx = zeros(Complex{T}, n_max + 1)
    cv = zeros(Complex{T}, n_max)
    dv = zeros(Complex{T}, n_max)
    logcv = zeros(Complex{T}, n_max)
    logdv = zeros(Complex{T}, n_max)
    logconjcv = zeros(Complex{T}, n_max)
    logconjdv = zeros(Complex{T}, n_max)
    cv_an = zeros(Complex{T}, n_max)
    dv_an = zeros(Complex{T}, n_max)
    

    # downward recursion
    @inbounds for n = n_max:-1:1
        rn = n + 1
        dx[n] = (rn / x) - (1 / (dx[rn] + rn / x))
    end
    
    gy0 = im #riccatibesselksi1(-1, y)/riccatibesselksi1(0, y)

    # upward recursion
    @inbounds for n = 1:nmax
        if n == 1 
            gy[n] = 1/(- gy0 + n / y ) - n/y
        else
            rn = n - 1 
            gy[n] = 1/(- gy[rn] + n / y ) - n/y
        end

        psi_n = riccatibesselpsi(n, x)
        ksi_n = riccatibesselksi1(n, x)
        psi_n_1 = riccatibesselpsi(n - 1, x)
        ksi_n_1 = riccatibesselksi1(n - 1, x)
        
       
        t_c = gy[n] + ref_idx * n / x
        t_d = gy[n] * ref_idx + n / x
        
        logcv[n] = log(-( t_c * ksi_n - ref_idx * ksi_n_1 )) - log( t_c * psi_n - ref_idx * psi_n_1 ) #TE
        logdv[n] = log(-( t_d * ksi_n - ksi_n_1 )) - log( t_d * psi_n - psi_n_1 ) # TM
        logconjcv[n] = log(conj(-( t_c * ksi_n - ref_idx * ksi_n_1 ))) - log(conj( t_c * psi_n - ref_idx * psi_n_1 )) #TE
        logconjdv[n] = log(conj(-( t_d * ksi_n - ksi_n_1 ))) - log(conj( t_d * psi_n - psi_n_1 )) # TM
    end
    return  logcv, logdv, logconjcv, logconjdv
end





struct SphericalShell{U,V} <: AbstractObject
    sphere::U
    cavity::V
end

sphere(spsh :: SphericalShell) = spsh.sphere
cavity(spsh :: SphericalShell) = spsh.cavity


function transmission_n(miecoeff_sphere, miecoeff_cavity, order)
    sph = miecoeff_sphere
    cav = miecoeff_cavity
    return (2 * order + 1) * (real(cav) + 1) * (real(sph) - abs(sph)^2) / abs(1 + cav * sph)^2 
end

function transmission_n(sphericalshell :: SphericalShell, w, order)
    (;sphere, cavity) = sphericalshell
    sph_te, sph_tm  = mie_coefficients(sphere, w, order)
    cav_te, cav_tm = mie_coefficients(cavity, w, order)
    tr_te = transmission_n(sph_te[order], cav_te[order], order)
    tr_tm = transmission_n(sph_tm[order], cav_tm[order], order)
    tr_tot = tr_te + tr_tm
    isnan(tr_tot) ? (@show w) : nothing
    return (te = tr_te, tm = tr_tm, total = tr_tot)  
end



function transmission(sphericalshell :: SphericalShell, w; nmax=0, rtol=1e-6)
    (;sphere, cavity) = sphericalshell
    xsph = size_parameter(sphere, w)
    xsph_stop, ysph_stop, nsph_max = max_order(sphere, w)
    
    xcav = size_parameter(cavity, w)
    xcav_stop, ycav_stop, ncav_max = max_order(cavity, w)
     
    scale = xcav-xsph

   n_max = min(nsph_max, ncav_max) #Integer(round(1/xcav/2)) + nsph_max
   n_max = Integer(round(xsph/(xcav - xsph)/3)) + n_max
    nmax != 0 ? n_max = nmax : nothing

    te_sph, tm_sph = mie_coefficients(sphere, w, n_max)
    te_cav, tm_cav = mie_coefficients(cavity, w, n_max)

    logte_sph, logtm_sph, logconjte_sph, logconjtm_sph = log_mie_coefficients(sphere, w, n_max)
    logte_cav, logtm_cav, logconjte_cav, logconjtm_cav = log_mie_coefficients(cavity, w, n_max)
    
    qw_te = 0.0
    qw_tm = 0.0
    
    qw_te_1 = 0.0
    qw_tm_1 = 0.0
    err = 100

    if nmax != 0
        for i in 1:nmax
            #qw_te += (2*i+1) * (real(te_cav[i]) + 1) * (real(te_sph[i]) - abs(te_sph[i])^2) / abs(1 + te_cav[i]*te_sph[i])^2
            #qw_tm += (2*i+1) * (real(tm_cav[i]) + 1) * (real(tm_sph[i]) - abs(tm_sph[i])^2) / abs(1 + tm_cav[i]*tm_sph[i])^2
            
            qw_te += (2*i+1) *(0.25 * exp(logte_cav[i] + logte_sph[i])
                    + 0.25 * exp(logconjte_cav[i] + logte_sph[i])
                    + 0.25 * exp(logte_cav[i] + logconjte_sph[i])
                    + 0.25 * exp(logconjte_cav[i] + logconjte_sph[i])
                    - 0.5 * exp(logte_cav[i] + logte_sph[i] + logconjte_sph[i])
                    - 0.5 * exp(logconjte_cav[i] + logte_sph[i] + logconjte_sph[i])
                    + real(te_sph[i]) 
                    - exp(logte_sph[i] + logconjte_sph[i]) )/ abs(1 + exp(logte_cav[i] + logte_sph[i]))^2
                    
            qw_tm += (2*i+1) *(0.25 * exp(logtm_cav[i] + logtm_sph[i])
                    + 0.25 * exp(logconjtm_cav[i] + logtm_sph[i])
                    + 0.25 * exp(logtm_cav[i] + logconjtm_sph[i])
                    + 0.25 * exp(logconjtm_cav[i] + logconjtm_sph[i])
                    - 0.5 * exp(logtm_cav[i] + logtm_sph[i] + logconjtm_sph[i])
                    - 0.5 * exp(logconjtm_cav[i] + logtm_sph[i] + logconjtm_sph[i])
                    + real(tm_sph[i]) 
                    - exp(logtm_sph[i] + logconjtm_sph[i]) )/ abs(1 + exp(logtm_cav[i] + logtm_sph[i]))^2
        
                end
        return (te =  real(qw_te), tm = real(qw_tm), total = real(qw_te + qw_tm), nmax = n_max)
    elseif nmax == 0
        i=0
        while (err > rtol) && (i < n_max)
            i += 1 
            #qw_te += (2*i+1) * (real(te_cav[i]) + 1) * (real(te_sph[i]) - abs(te_sph[i])^2) / abs(1 + te_cav[i]*te_sph[i])^2
            #qw_tm += (2*i+1) * (real(tm_cav[i]) + 1) * (real(tm_sph[i]) - abs(tm_sph[i])^2) / abs(1 + tm_cav[i]*tm_sph[i])^2       
    
            qw_te += real((2*i+1) *(0.25 * exp(logte_cav[i] + logte_sph[i])
                    + 0.25 * exp(logconjte_cav[i] + logte_sph[i])
                    + 0.25 * exp(logte_cav[i] + logconjte_sph[i])
                    + 0.25 * exp(logconjte_cav[i] + logconjte_sph[i])
                    - 0.5 * exp(logte_cav[i] + logte_sph[i] + logconjte_sph[i])
                    - 0.5 * exp(logconjte_cav[i] + logte_sph[i] + logconjte_sph[i])
                    + real(te_sph[i]) 
                    - exp(logte_sph[i] + logconjte_sph[i]) )/ abs(1 + exp(logte_cav[i] + logte_sph[i]))^2)
                    
            qw_tm += real((2*i+1) *(0.25 * exp(logtm_cav[i] + logtm_sph[i])
                    + 0.25 * exp(logconjtm_cav[i] + logtm_sph[i])
                    + 0.25 * exp(logtm_cav[i] + logconjtm_sph[i])
                    + 0.25 * exp(logconjtm_cav[i] + logconjtm_sph[i])
                    - 0.5 * exp(logtm_cav[i] + logtm_sph[i] + logconjtm_sph[i])
                    - 0.5 * exp(logconjtm_cav[i] + logtm_sph[i] + logconjtm_sph[i])
                    + real(tm_sph[i]) 
                    - exp(logtm_sph[i] + logconjtm_sph[i]) )/ abs(1 + exp(logtm_cav[i] + logtm_sph[i]))^2)
       
            
            err = abs(qw_te_1 + qw_tm_1 - qw_te - qw_tm)/abs(qw_te + qw_tm)
            qw_te_1 = qw_te
            qw_tm_1 = qw_tm
           
        end
        i > n_max ? @warn("The calculation did not converge. The number of iterations exceeds the value of n_max. Try increasing rtol.") : nothing 
        err > rtol ? @warn("The calculation did not converge. The relative error is larger than rtol. Try increasing rtol.") : nothing
        return (te =  qw_te, tm = qw_tm, total = (qw_te + qw_tm), numiter = i, nmax = n_max, rel_error = err)
    end
end


"Spectral heat flux"
function heat_flux(sphericalshell :: SphericalShell, w; nmax=0, rtol=1e-6)
    (;sphere, cavity) = sphericalshell
    T1 = temperature(sphere)
    T2 = temperature(cavity)
    delTheta = (bose_einstein(w,T1) - bose_einstein(w,T2))
    transm = transmission(sphericalshell, w; nmax, rtol)
    qw_te = delTheta * transm[:te] * 2 / pi
    qw_tm = delTheta * transm[:tm] * 2 / pi
    qw_tot = delTheta * transm[:total] * 2 / pi
    isnan(qw_tot) ? (@show w) : nothing
    return (te = qw_te, tm = qw_tm, total = qw_tot, nmax = transm[:nmax])
end

"Spectral heat flux for a single Mie order"
function heat_flux_n(sphericalshell :: SphericalShell, w, order ; rtol=1e-6)
    (;sphere, cavity) = sphericalshell
    T1 = temperature(sphere)
    T2 = temperature(cavity)
    delTheta = (bose_einstein(w,T1) - bose_einstein(w,T2))
    mie_sp_te, mie_sp_tm  = mie_coefficients(sphere, w, order)
    mie_cav_te, mie_cav_tm = mie_coefficients(cavity, w, order)
    transm_te = transmission_n(mie_sp_te[order], mie_cav_te[order], order)
    transm_tm = transmission_n(mie_sp_tm[order], mie_cav_tm[order], order)
    qw_te = delTheta * transm_te * 2 / pi
    qw_tm = delTheta * transm_tm * 2 / pi
    qw_tot = qw_te + qw_tm 
    isnan(qw_tot) ? (@show w) : nothing
    return (te = qw_te, tm = qw_tm, total = qw_tot)
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
    return (te = qte * kb / ħ, tm = qtm * kb / ħ , total = (qte + qtm) * kb / ħ, nmax = n_max_w)
end

"Total heat flux by mie order"
function heat_flux_n(sphericalshell :: SphericalShell, order; w1 = 0.0, w2 = Inf, rtol=1e-6, rtol_int = 1e-6)
    (;sphere, cavity) = sphericalshell
    hf_w(w) = heat_flux_n(sphericalshell, w, order; rtol)
    qte_w(u) = hf_w(u*kb/ħ)[:te]
    qtm_w(u) = hf_w(u*kb/ħ)[:tm]
    
    u1 = w1*ħ/kb
    u2 = w2*ħ/kb

    (qte, err) = quadgk(qte_w, u1 , u2 ; rtol = rtol_int)
    (qtm, err) = quadgk(qtm_w, u1 , u2 ; rtol = rtol_int)
    return (te = qte * kb / ħ, tm = qtm * kb / ħ , total = (qte + qtm) * kb / ħ)
end