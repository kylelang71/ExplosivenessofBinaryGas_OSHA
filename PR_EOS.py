import numpy as np


def molar_vapor_composition_STP(omega1, Tc1, Pc1, x1, omega2, Tc2, Pc2, x2):
    # Define States (OSHA LEL) and Constants
    T = 293.15  # K
    Tr1 = T / Tc1
    Tr2 = T / Tc2
    R = 8.314472  # Gas Constant in J / Mol*K
    P = 0.1013  # MPa

    # Define equations
    m1 = 0.37464 + 1.54226 * omega1 - 0.26992 * (omega1 ** 2)
    alpha1 = (1 + m1 * (1 - (Tr1 ** 0.5))) ** 2
    a11 = (0.45723553 * (R ** 2) * (Tc1 ** 2) / Pc1) * alpha1
    m2 = 0.37464 + 1.54226 * omega2 - 0.26992 * (omega2 ** 2)
    alpha2 = (1 + m2 * (1 - (Tr2 ** 0.5))) ** 2
    a22 = (0.45723553 * (R ** 2) * (Tc2 ** 2) / Pc2) * alpha2
    a12 = ((a11 * a22) ** 0.5) * (1 - 0.67704)
    """ Kij = 0.67704 provides 3.2% Err at 313.15K
    ["Correlation of Vapor-Liquid Equilibrium Data for Binary Mixtures
    Containing One or More Polar Component by the Parameters from Group
    Contributions Equation of State" - M Moshfegian & R. N. Maddox]"""
    a = (x1 ** 2 * a11) + 2 * (x1 * x2 * a12) + (x2 ** 2 * a22)
    coef = 0.07779607
    b1 = coef * R * Tc1 / Pc1
    b2 = coef * R * Tc2 / Pc2
    b = b1 * x1 + b2 * x2
    A = a * P / ((R ** 2) * (T ** 2))
    B = b * P / (R * T)

    # CREATE & SOLVE CUBIC (Elliot & Lira "Introductory Chemical Engineering Thermodynamics")
    M3 = 1
    M2 = -(1 - B)
    M1 = (A - 3 * (B ** 2) - 2 * B)
    M0 = -(A * B - B ** 2 - B ** 3)
    fxn = np.roots([M3, M2, M1, M0])

    # Calculate Component - Phase Fugacity Coefficients
    precalc_fcl1 = (b1 / b) * (fxn[2] - 1) - np.log(fxn[2] - B) - \
                   (A / (2 * (2 ** 0.5) * B) * (2 * (a11 / a) - (b1 / b))) * \
                   np.log(fxn[2] + (2.414 * B) / (fxn[2] - 0.414 * B))
    fc_liquid_1 = np.exp(precalc_fcl1)

    precalcfcv1 = (b1 / b) * (fxn[0] - 1) - np.log(fxn[0] - B) - \
                  (A / (2 * (2 ** 0.5) * B) * (2 * (a11 / a) - (b1 / b))) * \
                  np.log(fxn[0] + (2.414 * B) / (fxn[0] - 0.414 * B))
    fc_vapor_1 = np.exp(precalcfcv1)
    y1_calc = x1 * fc_liquid_1 / fc_vapor_1

    precalc_fcl2 = (b2 / b) * (fxn[2] - 1) - np.log(fxn[2] - B) - \
                   (A / (2 * (2 ** 0.5) * B) * (2 * (a22 / a) - (b2 / b))) * \
                   np.log(fxn[2] + (2.414 * B) / (fxn[2] - 0.414 * B))
    fc_liquid_2 = np.exp(precalc_fcl2)

    precalcfcv2 = (b2 / b) * (fxn[0] - 1) - np.log(fxn[0] - B) - \
                  (A / (2 * (2 ** 0.5) * B) * (2 * (a22 / a) - (b2 / b))) * \
                  np.log(fxn[0] + (2.414 * B) / (fxn[0] - 0.414 * B))
    fc_vapor_2 = np.exp(precalcfcv2)
    y2_calc = x1 * fc_liquid_2 / fc_vapor_2

    y = y2_calc + y1_calc
    y1 = y1_calc / y
    y2 = y2_calc / y

    return {'Comp1_MoleFraction': y1, 'Comp2_MoleFraction': y2}
