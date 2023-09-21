# Definitions for Mie coeffcients of cylinders (layered or not)
# Single cylinder, normal incidence is the implementation bu Bohren and Huffman

struct Cylinder{T <: OptProp, U <: Number } <: AbstractObject
    material :: T
    radius   :: U
end

#
function log_deriv!(g,y,n_max)
    for n = 1:(n_max-1)
        rn   = n_max - n + 1
        k    = n_max - n
        g[k] = (rn - 2)/y -  1/( g[k+1] + (rn - 1)/y )
     end
     return g
end

function _bhbesselj!(f,size_par,mst)
    f[mst + 1] = 0.0
    f[mst]     = 1e-32

    for l=1:(mst - 1)
        ml = mst - l
        f[ml] = 2.0*ml*f[ml+1]/size_par - f[ml+2]
    end
    # using F to compute alpha
    alpha = f[1]
    for ll=2:2:(mst-2)
        alpha += 2*f[ll + 1]
    end
    return f, alpha
end

function bhbesselj!(bj,f,size_par,mst)

    f,alpha = _bhbesselj!(f,size_par,mst)

    # From alpha and F, we compute J
    for n = 1:(mst - 1)
        bj[n] = f[n]/alpha
    end
    return bj
end

# We compute Y
function bhbessely!(by,bj,size_par,mst,x_stop)
    by[1] = bj[1]*(log(size_par/2) + γ)
    m = Int(mst/2) - 1
    for l=1:m
        by[1] = by[1] - 2*((-1)^l)*bj[2*l + 1]/l
    end
    by[1] = 2/pi * by[1]
    by[2] = bj[2]*by[1] - 2/pi/size_par
    by[2] = by[2]/bj[1]

    ns = round(Int,x_stop) - 1

    for kk = 1:ns
       by[kk+2] = 2*kk*by[kk+1]/size_par - by[kk]
    end
    return by
end


function mie_coefficient(cylinder :: Cylinder{T,U}, w) where {T<: OptProp, U<: Number}
    x_stop,y_stop,n_max = max_order(cylinder,w)
    size_par = size_parameter(cylinder,w)
    ref_idx  = refractive_index(cylinder,w)
    y        = size_par*ref_idx

    n_delta   = (101 + size_par)^(0.499)
    mst_temp  = round(Int,x_stop + n_delta)
    mst       = iseven(mst_temp) ? mst_temp : (mst_temp + 1)

    g  = zeros(Complex{U}, n_max)
    f  = zeros(U, mst+1)
    an = zeros(Complex{U}, round(Int,x_stop)+1)
    bn = zeros(Complex{U}, round(Int,x_stop)+1)

    bj = zeros(U, mst+1)
    by = zeros(U, round(Int,x_stop)+1)
    bh = zeros(Complex{U}, mst+1)

    # Computes G function for argument y
    g  = log_deriv!(g,y,n_max)

    # We compute J
    bj = bhbesselj!(bj,f,size_par,mst)

    # We compute Y
    by = bhbessely!(by,bj,size_par,mst,x_stop)

    # From J and Y, we obtain H
    nn = round(Int,x_stop) + 1
    for n=1:nn
        bh[n] = bj[n] + im*by[n]
    end

    # All Bessel function are computed, we can now compute de Mie coefficients
    a0 = g[1]*bj[1]/ref_idx + bj[2]
    a0 = a0/(g[1]*bh[1]/ref_idx + bh[2])
    b0 = ref_idx*g[1]*bj[1] + bj[2]
    b0 = b0/(ref_idx*g[1]*bh[1] + bh[2])
    qspar = abs(b0)^2
    qsper = abs(a0)^2

    for n=1:round(Int,x_stop)
        an[n] = (g[n+1]/ref_idx + n/size_par)*bj[n+1] - bj[n]
        an[n] = an[n]/((g[n+1]/ref_idx + n/size_par)*bh[n+1] - bh[n])

        bn[n] = (ref_idx*g[n+1] + n/size_par)*bj[n+1] - bj[n]
        bn[n] = bn[n]/((ref_idx*g[n+1] + n/size_par)*bh[n+1]-bh[n])
    end
    return [a0 ; an] , [b0 ; bn]
end


function efficiencies(cylinder :: Cylinder, w)
    size_par = size_parameter(cylinder,w)
    a,b      = mie_coefficient(cylinder, w)

    q_sca_par = abs(b[1])^2
    q_sca_per = abs(a[1])^2
    q_ext_par = real(b[1])
    q_ext_per = real(a[1])

    for n = 2:length(a)
        q_sca_par += 2*abs(b[n])^2
        q_sca_per += 2*abs(a[n])^2
        q_ext_par += 2*real(b[n])
        q_ext_per += 2*real(a[n])
    end

    return  2/(size_par)*q_sca_par,2/(size_par)*q_ext_par,2/(size_par)*q_sca_per,2/(size_par)*q_ext_per
end
