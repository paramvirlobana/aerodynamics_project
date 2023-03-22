import math

# Aircraft parameters
cbar_w = 1.66255  # Mean chord of the main wing (m)
s_w = 16.7225  # Wing area (m^2)
s_h = 3.34451  # Horizontal stabilizer area (m^2)
s_v = 1.105546  # Vertical stabilizer area (m^2)
b_w = 10.0584  # Span of the main wing (m)
b_h = 3.6576  # Span of the horizontal stabilizer (m)
h_v = 0.9144  # Distance above center of gravity to the aerodynamic center of the tail
c_lw_alpha = 4.44  # Lift slope of the main wing
c_lh_alpha = 3.97  # Lift slope of the horizontal stabilizer
c_lv_alpha = 3.40  # Lift slope of the horizontal stabilizer
l_w = -0.216408  # Distance aft of the center of gravity to aerodynamic center of the main wing (m)
l_h = 4.355592  # Distance aft of the center of gravity to aerodynamic center of the horizontal stabilizer (m)
l_v = 4.514088  # Distance aft of the center of gravity to aerodynamic center of the vertical stabilizer (m)
eta_h = 1.0  # Dynamic pressure ratio relative to the free stream on the horizontal stabilizer
eta_v = 1.0  # Dynamic pressure ratio relative to the free stream on the vertical stabilizer
eps_d_alpha = 0.44  # Down wash gradient, or the change in down wash with angle of attack
eps_s_beta_v = -0.10  # Side wash gradient, or the change in side wash with angle of attack
gamma = -0.1  # Wing dihedral angle (degrees)
kappa_l = 1.07  # Wing dihedral factor (See Figure 5.6.3 of Phillips)
kappa_gamma = 0.83  # Wing dihedral factor (See Figure 5.6.3 of Phillips)


def static_margin(s_w, s_h, c_lw_alpha, c_lh_alpha, l_w, l_h, eta_h, eps_d_alpha, cbar_w):
    # Example 4.4.1 of Phillips
    sm = (l_w * c_lw_alpha + (s_h * l_h) / (s_w) * eta_h * c_lh_alpha * (1 - eps_d_alpha)) / (
            cbar_w * (c_lw_alpha + s_h / s_w * eta_h * c_lh_alpha * (1 - eps_d_alpha)))

    # Print the static margin
    print('Static Margin: ', round(sm * 100, 2), '% (recommended 5->15%)', sep='')

    return


def yaw_derivative(eta_v, s_v, s_w, l_v, l_w, b_w, c_lv_alpha, eps_s_beta_v):
    # Equation 5.2.7 of Phillips
    deltac_n_beta_v = eta_v * (s_v * l_v) / (s_w * b_w) * c_lv_alpha * (1 - eps_s_beta_v)

    # Print the yaw stability derivative
    print('Yaw Stability Derivative: ', round(deltac_n_beta_v, 3),
          ' (recommended 0.06->0.15, vertical stabilizer contribution only!)', sep='')

    return


def roll_derivative(gamma, kappa_gamma, kappa_l, c_lw_alpha, h_v):
    # Convert the dihedral angle to radians
    gamma = gamma * math.pi / 180

    # Get the contribution from the main wing, assuming negligible sweep, via Equation 5.6.13 in Phillips
    deltac_l_beta_gammaw = -(2 * math.sin(gamma)) / (
            3 * math.pi * math.cos(gamma) ** 4) * kappa_gamma * kappa_l * c_lw_alpha

    # Get the contribution from the vertical stabilizer via Equation 5.6.22 in Phillips
    deltac_l_beta_v = -eta_v * (s_v * h_v) / (s_w * b_w) * (1 - eps_s_beta_v) * c_lv_alpha

    # Get the contribution from the horizontal stabilizer via Equation 5.6.23 in Phillips (assuming conventional tail)
    deltac_l_beta_h = + 0.08 * eta_v * (s_v * b_h) / (s_w * b_w) * (1 - eps_s_beta_v) * c_lv_alpha

    # Get the total roll stability derivative
    deltac_n_beta_v = deltac_l_beta_gammaw + deltac_l_beta_v + deltac_l_beta_h

    # Print the roll stability derivative
    print('Roll Stability Derivative: ', round(deltac_n_beta_v, 3), ' (recommended -0.1->0)', sep='')

    return


# Compute the static margin for longitudinal stability
static_margin(s_w, s_h, c_lw_alpha, c_lh_alpha, l_w, l_h, eta_h, eps_d_alpha, cbar_w)

# Compute the yaw stability derivative
yaw_derivative(eta_v, s_v, s_w, l_v, l_w, b_w, c_lv_alpha, eps_s_beta_v)

# Compute the roll stability derivative
roll_derivative(gamma, kappa_gamma, kappa_l, c_lw_alpha, h_v)
