
module LPC

#
# TODO:
# [v] lpc-burg
# [ ] lattice-filter synth
# [ ] warped-lpc
# [ ] lpc-acorr
# [ ] lpc-covar
# [ ] Line-Spectral-Pairs
# [v] power-spectrum est
#

export lpc_burg, pow_spect
export log_area_ratio, inv_sine_coeffs

function lpc_burg{T <: Number}( x::AbstractVector{T}, p::Int)
# LPC (Linear-Predictive-Code) estimation, using the Burg-method
# @param: x - signal to predict
# @param: p - prediction order
# @retval: a - direct form coeffs
# @retval: r - reflection coefficients (usable for lattice-filter?)
# such that:
# minimizing x_prediction(n) = a*x
# ie: an IIR:
#  y = filt(1, a, white)
# should have similiar spectrum to x
# where 'white' is a signal with white spectrum
# 
# The Burg-Method, which, provide both a minimal-phase (stable IIR)
# and estimates 'a', on a finite support [1]
# Combining the advantages on both traditional LPC methods
# (autocorrelation - minPhase, and covariance - finite-support)
#
# Implementing the loop from the article:
# [1] - ENHANCED PARTIAL TRACKING USING LINEAR PREDICTION
# (DAFX'03 article, Lagrange et al)
# http://www.elec.qmul.ac.uk/dafx03/proceedings/pdfs/dafx19.pdf
#
# Note: After implementing - I found matlab has an 'arburg()' - in
# the 'signal-processing-toolbox', with virtually identical code
# but hey, why use a library when you can easily reimplement it?
# anyhow, normalized interface.
#
# (c) Kobic, 2009, 2015, MIT License
# Translated my old(er) matlab code into julia::

# TODO: check memory allocations, gc etc!

    ef = x    # forward error
    eb = x    # backwards error
    a = [1; zeros(T, p)]  # prediction coeffs
    refl = zeros(T, p) # reflection coeffs
    k = zero(T)

# zero-mean - or not?
    prediction_err = dot(x,x) ./ length(x)  # variance

    for m in 1:p
        efp = ef[1:end-1]
        ebp = eb[2:end]
        k = -2 .* dot(ebp,efp) ./ ( dot(ebp,ebp) .+ dot(efp,efp) )
        refl[m] = k           # save reflection coeffs (if needed)
        ef = efp + k.*ebp
        eb = ebp + k.*efp
        a[1:m+1] = [a[1:m]; 0] + k.*[0; a[m:-1:1]]
        prediction_err *= (1 - k*k)
    end

    return a, prediction_err, refl 
end



function pow_spect(ar, err, dft_order = 256)
    f = rfft( [ar; zeros(dft_order - length(ar))] )
    psd = err ./ (abs(f).^2 + eps() )

    return psd
end

function log_area_ratio(k)
    return log( (1 + k) ./ (1 - k) )
end

function inv_sine_coeffs(k)
    return 2asin(k) ./ pi
end

end # module

