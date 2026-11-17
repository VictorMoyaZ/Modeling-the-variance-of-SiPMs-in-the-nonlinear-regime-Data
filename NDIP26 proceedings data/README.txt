DATA FOR "AN ANALYTICAL CHARGE DISPERSION CORRECTION TO NON-PARALYZABLE
DEAD-TIME MODELS FOR PASSIVE SIPM VARIANCE"
================================================================================

Data.xlsx has one sheet per pulse duration: T_pulse/tau_rec = 0.01, 1, 10 and
100 (tau_rec = 30 ns). Equation numbers refer to the manuscript.


COLUMNS
--------------------------------------------------------------------------------

Simulated quantities:

  x              Expected number of photoelectrons per pixel (simulation input).
  y              E[Q]/(N*q): mean output charge per pixel, in units of the
                 avalanche charge.
  Std y          Standard deviation of y over the Monte Carlo repetitions.
                 The standard error of the mean is (Std y)/sqrt(n_rep).

Parameters fitted by least squares (given in the first row of each sheet only):

  Tdead/Tpulse   Dead-time model, Eq. (4): y = x/(1 + x*Tdead/Tpulse).
                 Also the fitted parameter of Eq. (15).
  a, b           Partial pixel recovery model, Eq. (9) with c = d = 0.
  alpha, beta    Charge-resolution model, Eq. (11).

Other parameters (first row only):

  epsilon        Parameter of Eq. (11). Equal to 1 in all sheets.
  1+(s/q)^2      1 + s^2/q^2 = 1.0225, with s/q = 0.15 (fixed, not fitted).


MONTE CARLO REPETITIONS (n_rep)
--------------------------------------------------------------------------------

  x <= 1            10000
  1 < x <= 10        5000
  x > 10            1000
