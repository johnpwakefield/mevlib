#!/usr/bin/env julia


using LinearAlgebra


k0s = [         # 1 / s
    0.000 1.413 4.337 1.163 0.114 0.386;
    0.000 0.000 0.229 0.161 0.041 0.137;
    0.000 0.000 0.000 0.128 0.030 0.103;
    0.000 0.000 0.000 0.000 0.000 0.000;
    0.000 0.000 0.000 0.000 0.000 0.000;
    0.000 0.000 0.000 0.000 0.000 0.000
       ]

Eas = [         # kJ / mol
    00.0 47.6 43.4 38.5 30.2 30.0;
    00.0 00.0 54.1 62.9 66.7 65.0;
    00.0 00.0 00.0 80.5 85.2 77.3;
    00.0 00.0 00.0 00.0 00.0 00.0;
    00.0 00.0 00.0 00.0 00.0 00.0;
    00.0 00.0 00.0 00.0 00.0 00.0
       ]

Ms = [          # g / mol
      444.0, 230.0, 115.0, 52.0, 16.0, 400.0
      ]


T  = 673                    # Kelvin
T0 = 773                    # Kelvin
R = 8.31446261815324e-3     # kJ per mol Kelvin
Dpore = 2.0                 # nm
epspore = 0.319             # nondim
tau = 7.0                   # nondim


# the units in this formula do work out correctly
Dis = @. 48.5 * Dpore * sqrt(T / Ms) * epspore / tau

function kij(k0, Ea, R, T, T0)
    return k0 * exp(-Ea / R * (1 / T - 1 / T0))
end


kijs = kij.(k0s, Eas, R, T, T0)


nB = -transpose((Diagonal([sum(k0s[i, :]) for i in 1:6]) - k0s) ./ Dis)
R, Lam = eigvecs(nB), diag(nB)


bdry = [1.0, 0.8, 0.6, 0.4, 0.2, 0.0]
equilibrium = R * Diagonal([Lam[i] == 0.0 ? 1.0 : 0.0 for i in 1:6]) * (R \ bdry)


println("\nk0s")
display(transpose(k0s))
println("\nEas")
display(transpose(Eas))
println("\nDis")
display(Dis)
println("\nkijs")
display(kijs)
println("\nnB")
display(nB)
println("\nR")
display(R)
println("\nequilibrium")
display(equilibrium)
println("\ninitsum = $(sum(bdry)), eqsum = $(sum(equilibrium))")
println("\n")


