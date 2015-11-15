# LPC.jl
Linear predictive coding, Julia module

a, e, rc = lpc_burg( x, p)

x - input signal, p - prediction order

a - predicted all-pole filter
e - variance of prediction error
rc - reflection-coefficients 

Currently implemented the Burg method to calculate lpc, which calculates the Reflection-Coefficients as a side-effect, and has several stability advantedges.


 LPC (Linear-Predictive-Code) estimation, using the Burg-method
 @param: x - signal to predict
 @param: p - prediction order
 @retval: a - direct form coeffs
 @retval: r - reflection coefficients (usable for lattice-filter?)
 such that:
 minimizing x_prediction(n) = a*x
 ie: an IIR:
  y = filt(1, a, white)
 should have similiar spectrum to x
 where 'white' is a signal with white spectrum
 
 The Burg-Method, which, provide both a minimal-phase (stable IIR)
 and estimates 'a', on a finite support [1]
 Combining the advantages on both traditional LPC methods
 (autocorrelation - minPhase, and covariance - finite-support)

 Implementing the loop from the article:
 [1] - ENHANCED PARTIAL TRACKING USING LINEAR PREDICTION
 (DAFX'03 article, Lagrange et al)
 http://www.eecs.qmul.ac.uk/legacy/dafx03/proceedings/pdfs/dafx19.pdf
