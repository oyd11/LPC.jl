# LPC.jl
Linear predictive coding, Julia module

lpc_burg( x, p)

x - input signal, p - prediction order

Currently implemented the Burg method to calculate lpc, which calculates the Reflection-Coefficients as a side-effect, and has several stability advantedges.

