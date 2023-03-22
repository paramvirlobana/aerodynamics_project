import numpy as np
import math

rho = 1.22  # Air density at current altitude (kg/m^3)
s_w = 17.1871  # Wing area (m^2)
b_w = 10.06  # Span of the main wing (m)
m = 1270.06  # Aircraft mass (kg)
g = 9.805  # Acceleration due to gravity (m/s^2)
v_0 = 54.86  # Airspeed (m/s)
c_d0 = 0.05  # Parasitic drag coefficient
i_xx_b = 1355.81  # Moment of inertia
i_yy_b = 4067.43  # Moment of inertia
i_zz_b = 4745.33  # Moment of inertia
i_xz_b = 40.67  # Moment of inertia
c_l_alpha = 4.40  # Lift slope of the aircraft
c_d_alpha = 0.35  # Drag slope of the aircraft
c_m_alpha = -0.68  # Moment slope of the aircraft
c_l_alphahat = 1.60
c_m_alphahat = -4.35
c_y_beta = -0.560
c_l_beta = -0.075
c_n_beta = 0.070
c_d_qbar = 0.0
c_l_qbar = 3.80
c_m_qbar = -9.95
c_y_pbar = 0.0
c_l_pbar = -0.410
c_n_pbar = -0.0575
c_y_rbar = 0.240
c_l_rbar = 0.105
c_n_rbar = -0.125
theta_0 = 0  # (rad)

# Pre-calculations
w = m * g


def longitudinal_stability(s_w, b_w, g, v_0, c_d0, i_yy_b, c_l_alpha, c_d_alpha, c_m_alpha, c_l_alphahat, c_m_alphahat,
                           c_l_qbar, c_m_qbar, w, theta_0):
    # Get the mean aerodynamic chord length
    cbar_w = s_w / b_w

    # Get the initial lift coefficient in steady level flight
    c_l0 = w * 1 / (0.5 * rho * v_0 ** 2 * s_w)

    # Get the matrix entries
    r_gx = g * cbar_w / (2 * v_0 ** 2)
    r_z_alphahat = - rho * s_w * cbar_w / (4 * w / g) * c_l_alphahat
    r_m_alphahat = rho * s_w * cbar_w ** 3 / (8 * i_yy_b) * c_m_alphahat
    r_x_mu = - rho * s_w * cbar_w / (4 * w / g) * (2 * c_d0)
    r_z_mu = - rho * s_w * cbar_w / (4 * w / g) * (2 * c_l0)
    r_m_mu = 0
    r_x_alpha = rho * s_w * cbar_w / (4 * w / g) * (c_l0 - c_d_alpha)
    r_z_alpha = rho * s_w * cbar_w / (4 * w / g) * (- c_l_alpha - c_d0)
    r_m_alpha = rho * s_w * cbar_w ** 3 / (8 * i_yy_b) * c_m_alpha
    r_x_qbar = 0.0
    r_z_qbar = - rho * s_w * cbar_w / (4 * w / g) * (c_l_qbar)
    r_m_qbar = rho * s_w * cbar_w ** 3 / (8 * i_yy_b) * c_m_qbar

    A1 = np.zeros((6, 6));
    A2 = np.zeros((6, 6));

    A1[0, 0] = r_x_mu
    A1[0, 1] = r_x_alpha
    A1[0, 2] = r_x_qbar
    A1[0, 5] = -r_gx * math.cos(theta_0)
    A1[1, 0] = r_z_mu
    A1[1, 1] = r_z_alpha
    A1[1, 2] = (1 + r_z_qbar)
    A1[1, 5] = - r_gx * math.sin(theta_0)
    A1[2, 0] = r_m_mu
    A1[2, 1] = r_m_alpha
    A1[2, 2] = r_m_qbar
    A1[3, 0] = math.cos(theta_0)
    A1[3, 1] = math.sin(theta_0)
    A1[3, 5] = - math.sin(theta_0)
    A1[4, 0] = - math.sin(theta_0)
    A1[4, 1] = math.cos(theta_0)
    A1[4, 5] = - math.cos(theta_0)
    A1[5, 2] = 1

    A2[0, 0] = 1
    A2[1, 1] = (1 - r_z_alphahat)
    A2[2, 1] = - r_m_alphahat
    A2[2, 2] = 1
    A2[3, 3] = 1
    A2[4, 4] = 1
    A2[5, 5] = 1

    # Solve the eigenvalue problem
    A = np.dot(np.linalg.inv(A2), A1)
    eigvals, eigvecs = np.linalg.eig(A)

    # Extract the largest and smallest eigenvalues
    eigmax = 0
    eigmin = 1E12
    for i in range(0, eigvals.shape[0]):
        if (np.abs(eigvals[i]) > np.abs(eigmax)):
            eigmax = eigvals[i]
        if ((np.abs(eigvals[i]) < np.abs(eigmin)) and (np.abs(eigvals[i]) > 1E-10)):
            eigmin = eigvals[i]

    # Calculate the short period mode
    sigma = - np.real(eigmax) * 2 * v_0 / cbar_w

    # Get the short period 99% damping time
    t_99 = - math.log(0.01) / sigma

    # Get the damped period
    omega_d = np.abs(np.imag(eigmax)) * 2 * v_0 / cbar_w
    period_d = 2 * math.pi / omega_d

    print('Damped Short-Period Time: ', round(period_d, 3), 's, 99% Damping: ', round(t_99, 3), 's', sep='')

    # Calculate the long period phugoid mode
    sigma = - np.real(eigmin) * 2 * v_0 / cbar_w

    # Get the long period 99% damping time
    t_99 = - math.log(0.01) / sigma

    # Get the damped period
    omega_d = np.abs(np.imag(eigmin)) * 2 * v_0 / cbar_w
    period_d = 2 * math.pi / omega_d

    print('Damped Long-Period (Phugoid) Time: ', round(period_d, 3), 's, 99% Damping: ', round(t_99, 3), 's', sep='')


def lateral_stability(s_w, b_w, g, v_0, i_xx_b, i_zz_b, i_xz_b, c_y_beta, c_l_beta, c_n_beta, c_y_pbar, c_l_pbar,
                      c_n_pbar, c_y_rbar, c_l_rbar, c_n_rbar, w, theta_0):
    r_y_beta = rho * s_w * b_w / (4 * w / g) * c_y_beta
    r_l_beta = rho * s_w * b_w ** 3 / (8 * i_xx_b) * c_l_beta
    r_n_beta = rho * s_w * b_w ** 3 / (8 * i_zz_b) * c_n_beta
    r_y_pbar = rho * s_w * b_w / (4 * w / g) * c_y_pbar
    r_l_pbar = rho * s_w * b_w ** 3 / (8 * i_xx_b) * c_l_pbar
    r_n_pbar = rho * s_w * b_w ** 3 / (8 * i_zz_b) * c_n_pbar
    r_y_rbar = rho * s_w * b_w / (4 * w / g) * c_y_rbar
    r_l_rbar = rho * s_w * b_w ** 3 / (8 * i_xx_b) * c_l_rbar
    r_n_rbar = rho * s_w * b_w ** 3 / (8 * i_zz_b) * c_n_rbar
    r_gy = g * b_w / (2 * v_0 ** 2)
    t_xz = i_xz_b / i_xx_b
    t_zx = i_xz_b / i_zz_b

    A1 = np.zeros((6, 6));
    A2 = np.zeros((6, 6));

    A1[0, 0] = r_y_beta
    A1[0, 1] = r_y_pbar
    A1[0, 2] = (r_y_rbar - 1)
    A1[0, 4] = r_gy * math.cos(theta_0)
    A1[1, 0] = r_l_beta
    A1[1, 1] = r_l_pbar
    A1[1, 2] = r_l_rbar
    A1[2, 0] = r_n_beta
    A1[2, 1] = r_n_pbar
    A1[2, 2] = r_n_rbar
    A1[3, 0] = 1
    A1[3, 5] = math.cos(theta_0)
    A1[4, 1] = 1
    A1[4, 2] = math.tan(theta_0)
    A1[5, 2] = 1 / math.cos(theta_0)

    A2[0, 0] = 1
    A2[1, 1] = 1
    A2[2, 2] = 1
    A2[3, 3] = 1
    A2[4, 4] = 1
    A2[5, 5] = 1
    A2[1, 2] = - t_xz
    A2[2, 1] = - t_zx

    # Solve the eigenvalue problem
    A = np.dot(np.linalg.inv(A2), A1)
    eigvals, eigvecs = np.linalg.eig(A)

    # Calculate the Dutch roll mode
    for i in range(0, eigvals.shape[0]):
        if np.abs(np.imag(eigvals[i])) > 1E-10:
            lam = eigvals[i]

    sigma = - np.real(lam) * 2 * v_0 / b_w
    t_99 = - math.log(0.01) / sigma

    omega_d = np.abs(np.imag(lam)) * 2 * v_0 / b_w
    period_d = 2 * math.pi / omega_d

    print('Damped Dutch Roll Time: ', round(period_d, 3), 's, 99% Damping: ', round(t_99, 3), 's', sep='')

    # Calculate the roll mode
    lam = 0
    for i in range(0, eigvals.shape[0]):
        if (np.abs(np.imag(eigvals[i])) < 1E-10 and np.abs(np.real(eigvals[i])) > np.abs(lam)):
            lam = np.real(eigvals[i])

    sigma = - np.real(lam) * 2 * v_0 / b_w
    t_99 = - math.log(0.01) / sigma

    print('Roll Mode 99% Damping: ', round(t_99, 3), 's', sep='')

    # Calculate the spiral mode
    lam = 1E12
    for i in range(0, eigvals.shape[0]):
        if (np.abs(np.imag(eigvals[i])) < 1E-10 and np.abs(np.real(eigvals[i])) < np.abs(lam) and np.abs(
                np.real(eigvals[i])) > 1E-10):
            lam = np.real(eigvals[i])

    sigma = - np.real(lam) * 2 * v_0 / b_w
    t_99 = - math.log(0.01) / sigma

    print('Spiral Mode 99% Damping: ', round(t_99, 3), 's', sep='')


longitudinal_stability(s_w, b_w, g, v_0, c_d0, i_yy_b, c_l_alpha, c_d_alpha, c_m_alpha, c_l_alphahat, c_m_alphahat,
                       c_l_qbar, c_m_qbar, w, theta_0)

lateral_stability(s_w, b_w, g, v_0, i_xx_b, i_zz_b, i_xz_b, c_y_beta, c_l_beta, c_n_beta, c_y_pbar, c_l_pbar, c_n_pbar,
                  c_y_rbar, c_l_rbar, c_n_rbar, w, theta_0)
