import numpy as np

# Gravitational parameter for Earth in km^3/s^2
GM = 3.986004418E+5

def calculate_angular_momentum(r_vector, v_vector):
    return np.cross(r_vector, v_vector)

def calculate_inclination(c_vector):
    c_z = c_vector[2]
    c = np.linalg.norm(c_vector)
    i = np.arccos(c_z/c)
    return i

def calculate_semi_major_axis(r_vector, v_vector):
    # Calculate semi-major axis
    v = np.linalg.norm(v_vector) # Norm of initial velocity
    r = np.linalg.norm(r_vector) # Norm of initial coordinates (radius)
    h = (v**2)/2 - GM/r # Orbital energy

    a = -GM / (2 * h)
    return a

def calculate_laplace_runge_lenz_vector(v_vector, c_vector, r_vector):
    r = np.linalg.norm(r_vector)
    l_vector = np.cross(v_vector, c_vector) - (GM/r) * r_vector
    
    return l_vector

def calculate_eccentricity(l_vector):
    l = np.linalg.norm(l_vector)
    return l / GM

def calculate_sin_cos_longitude_ascending_node(c_vector, i):
    c = np.linalg.norm(c_vector)
    # Calculate sin(Ω)
    c1 = c_vector[0]
    sin_Omega = c1 / (c * np.sin(i))
    
    # Calculate cos(Ω)
    c2 = c_vector[1]
    cos_Omega = - c2 / (c * np.sin(i))

    return (sin_Omega, cos_Omega)

def calculate_longitude_ascending_node(c_vector, i):
    sin_Omega, cos_Omega = calculate_sin_cos_longitude_ascending_node(c_vector, i)
    Omega = np.arctan2(sin_Omega, cos_Omega) # In Radians
    Omega_deg = np.degrees(Omega)
    
    return (Omega, Omega_deg)

def calculate_sin_cos_argument_pericenter(l_vector, c_vector, i):
    l1 = l_vector[0]
    l2 = l_vector[1]
    l3 = l_vector[2]
    l = np.linalg.norm(l_vector)

    # Calculate sin(ω)
    sin_omega = l3 / (l * np.sin(i))

    sin_Omega, cos_Omega = calculate_sin_cos_longitude_ascending_node(c_vector, i)

    # Calculate cos(ω)
    cos_omega = l1 / l * cos_Omega + l2 / l * sin_Omega
    
    return (sin_omega, cos_omega)

def calculate_argument_pericenter(l_vector, c_vector, i):
    sin_omega, cos_omega = calculate_sin_cos_argument_pericenter(l_vector, c_vector, i)
    
    omega = np.arctan2(sin_omega, cos_omega)
    omega_deg = np.degrees(omega)
    
    return (omega, omega_deg)
    
def calculate_argument_latitude(r_vector, c_vector, i):
    x = r_vector[0]
    y = r_vector[1]
    z = r_vector[2]
    r = np.linalg.norm(r_vector)
    
    sin_Omega, cos_Omega = calculate_sin_cos_longitude_ascending_node(c_vector, i)
    
    sin_u = z / (r * np.sin(i))
    cos_u = x / r * cos_Omega + y / r * sin_Omega
    
    return (sin_u, cos_u)

def calculate_sin_cos_true_anomaly(r_vector, l_vector, c_vector, i):
    sin_u, cos_u = calculate_argument_latitude(r_vector, c_vector, i)
    sin_omega, cos_omega = calculate_sin_cos_argument_pericenter(l_vector, c_vector, i)
    
    sin_nu = sin_u * cos_omega - cos_u * sin_omega
    cos_nu = cos_u * cos_omega + sin_u * sin_omega
    
    return sin_nu, cos_nu

def calculate_sin_eccetnric_anomaly_eccentric_anomaly(r_vector, l_vector, c_vector, i):
    e = calculate_eccentricity(l_vector)
    sin_nu, cos_nu = calculate_sin_cos_true_anomaly(r_vector, l_vector, c_vector, i)
    
    sin_E = (np.sqrt(1 - e**2) * sin_nu) / (1 + e * cos_nu)
    cos_E = (cos_nu + e) / (1 + e * cos_nu)
    
    E = np.arctan2(sin_E, cos_E)
    
    return (sin_E, E)

def calculate_mean_anomaly(r_vector, l_vector, c_vector, i):
    e = calculate_eccentricity(l_vector)
    sin_E, E = calculate_sin_eccetnric_anomaly_eccentric_anomaly(r_vector, l_vector, c_vector, i)
    
    M = E - e * sin_E
    M_deg = np.degrees(M)
    
    return M, M_deg

def orbital_params(initial_coordinates, initial_velocity):
    # Calculate necessary parameters
    c_vector = calculate_angular_momentum(initial_coordinates, initial_velocity)
    l_vector = calculate_laplace_runge_lenz_vector(initial_velocity, c_vector, initial_coordinates)

    # Calculate each orbital element using the provided functions
    semi_major_axis = calculate_semi_major_axis(initial_coordinates, initial_velocity)
    eccentricity = calculate_eccentricity(l_vector)
    inclination = calculate_inclination(c_vector)
    Omega_rads, Omega_degs = calculate_longitude_ascending_node(c_vector, inclination)
    omega_rads, omega_degs = calculate_argument_pericenter(l_vector, c_vector, inclination)
    mean_anomaly_rads, mean_anomaly_degs = calculate_mean_anomaly(initial_coordinates, l_vector, c_vector, inclination)
    
    return semi_major_axis, eccentricity, inclination, Omega_rads, Omega_degs, omega_rads, omega_degs, mean_anomaly_rads, mean_anomaly_degs

def critical_argument_1(u, m, M, Omega, omega, theta):
    return u * (M + omega + Omega) - m * theta

def critical_argument_2(u, m, M, Omega, omega, theta):
    return u * (M + omega) + m * (Omega - theta)

def critical_argument_3(u, m, M, Omega, omega, theta):
    return u * M + m * (omega + Omega - theta)

def critical_argument_4(u, m, M, Omega, omega, theta):
    return u * (M - Omega + omega) - m * theta

def critical_argument_5(u, m, M, Omega, omega, theta):
    return u * M + m * (-omega + 2 * Omega - theta)

def calculate_critical_arguments(u, m, M, Omega, omega, theta):
    phi1 = critical_argument_1(u, m, M, Omega, omega, theta)
    phi2 = critical_argument_2(u, m, M, Omega, omega, theta)
    phi3 = critical_argument_3(u, m, M, Omega, omega, theta)
    phi4 = critical_argument_4(u, m, M, Omega, omega, theta)
    phi5 = critical_argument_5(u, m, M, Omega, omega, theta)
    
    return phi1, phi2, phi3, phi4, phi5

# h in radians
def calculate_sidereal_time(jd):
    jd2000 = 2451545 # Julian date for J2000.0
    jdyear = 36525 # Days in Julian century

    # Fractional part of Julian day
    m = jd % 1 - 0.5

    # Days since J2000.0
    d = jd - m - jd2000

    # Time in Julian centuries since J2000.0
    t = (d + m) / jdyear

    # Fractional day in seconds
    mm = m * 86400

    # Sidereal time
    # d* = d + m
    # H0 = 244110.54841 + 8640184.812866 * (d* / 36525) + 0.093104 * t^2 - 6.2E-6 * t^3 
    s = (24110.54841 + mm + 236.555367908 * (d + m) + (0.093104 * t - 6.21e-6 * t**2) * t) # In seconds
    s_rad = s / 86400 * 2 * np.pi

    return s_rad
