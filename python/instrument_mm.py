from pyMuellerMat import common_mms as cmm
from pyMuellerMat import MuellerMat
import stokes
import angles
import numpy as np

def MODHIS_full_system_mm(pa = 0, altitude = 0, delta_HWP = 0.5, HWP_ang = 0, 
    wollaston_beam = 'o', TMT_matrix_noise = 0, NFIRAOS_matrix_noise = 0, MODHIS_matrix_noise = 0):
    """
    Returns the Mueller matrix of M3 with rotation.

    Args:
        delta_m3: (float) retardance of M3 (waves)
        epsilon_m3: (float) diattenuation of M3
        parang: (float) parallactic angle (degrees)
        altitude: (float) altitude angle in header (degrees)
        offset: (float) offset angle of M3 (degrees) - fit from M3 diattenuation fits
        TMT_matrix_noise: (float) % noise added to all Mueller matrix elements for TMT mirrors
        NFIRAOS_matrix_noise: (float) % noise added to all Mueller matrix elements for NFIRAOS
        MODHIS_matrix_noise: (float) % noise added to all Mueller matrix elements for MODHIS
    """

    # Parallactic angle rotation
    parang_rot = cmm.Rotator(name = "parang")
    parang_rot.properties['pa'] = pa

    # All telescope mirrors
    TMT = cmm.ArbitraryMatrix(name = "TMT")
    TMT.properties['mm'] = TMT_matrix(TMT_matrix_noise = TMT_matrix_noise)

     # TODO: Verify if this is positive or negative
    # Altitude angle rotation
    alt_rot = cmm.Rotator(name = "altitude")
    alt_rot.properties['pa'] = altitude

    # HWP for differencing - assumed to be before NFIRAOS
    hwp = cmm.Retarder(name = 'hwp') 
    hwp.properties['phi'] = delta_HWP * 2 * np.pi
    hwp.properties['theta'] = HWP_ang

    # NFIRAOS matrix
    NFIRAOS = cmm.ArbitraryMatrix(name = "NFIRAOS")
    NFIRAOS.properties['mm'] = NFIRAOS_matrix(NFIRAOS_matrix_noise = NFIRAOS_matrix_noise)

    # MODHIS matrix
    MODHIS = cmm.ArbitraryMatrix(name = "MODHIS")
    MODHIS.properties['mm'] = MODHIS_matrix(MODHIS_matrix_noise = MODHIS_matrix_noise)

    wollaston = cmm.WollastonPrism()
    wollaston.properties['beam'] = wollaston_beam
    
    sys_mm = MuellerMat.SystemMuellerMatrix([wollaston, MODHIS, parang_rot, 
        NFIRAOS, hwp, alt_rot, TMT, parang_rot])
    inst_matrix = sys_mm.evaluate()

    # # For debugging Wollaston + HWP matrix
    wollaston_and_HWP_mm = MuellerMat.SystemMuellerMatrix([wollaston, hwp])
    wollaston_and_HWP_matrix = wollaston_and_HWP_mm.evaluate()

    return inst_matrix, wollaston_and_HWP_matrix

# def full_system_mm_single_diff(
#     pa = 0, altitude = 0, delta_HWP = 0.5, HWP_ang = 0, TMT_matrix_noise = 0, 
#     NFIRAOS_matrix_noise = 0, MODHIS_matrix_noise = 0 factor = 1,
#     change_first_I_term = False):
#     """
#     Calculates an instrument matrix for the double difference or double
#     sum

#     NOTE: See Boris' overleaf file "VAMPIRES Integral Pol" for more details
    
#     Args:
#         fixed_params: (list) 

#     Returns:
#         data: (np.array) np.array([double_diff_matrix, double_sum_matrix])
#     """
#     # print("Fixed Params:  " + str(fixed_params))

#     FL1 = model(*fixed_params, parang, altitude, 
#                                      HWP_ang, IMR_ang, 1, 1)
#     FR1 = model(*fixed_params, parang, altitude, 
#                                      HWP_ang, IMR_ang, 2, 1)
#     FL2 = model(*fixed_params, parang, altitude, 
#                                      HWP_ang, IMR_ang,  1, 2)
#     FR2 = model(*fixed_params, parang, altitude, 
#                                      HWP_ang, IMR_ang,  2, 2)
    
#     # print("FL1: " + str(FL1))
#     # print("FR1: " + str(FR1))
#     # print("FL2: " + str(FL2))
#     # print("FR2: " + str(FR2))

#     double_diff_matrix = ((FL1 - FR1) - (FL2 - FR2)) / factor
#     double_sum_matrix = ((FL1 + FR1) + (FL2 + FR2)) / factor

#     if change_first_I_term:
#         double_diff_matrix[0, 0] = 1

#     return np.array([double_diff_matrix, double_sum_matrix])

def vary_full_system_mm(mueller_matrix, noise_percentage=0):
    """
    Adds Gaussian noise to each element of a Mueller matrix.

    Parameters:
    ----------
    mueller_matrix : np.ndarray
        The input Mueller matrix (assumed to be a 4x4 matrix).
    noise_percentage : float, optional
        The percentage of Gaussian noise to add to each element (default is 1%).

    Returns:
    -------
    varied_mueller_matrix : np.ndarray
        The Mueller matrix with added noise.
    """
    # Calculate the noise as a fraction of each element
    noise = np.random.normal(0, noise_percentage / 100, mueller_matrix.shape)
    varied_mueller_matrix = mueller_matrix * (1 + noise)
    
    return varied_mueller_matrix

def propagate_on_sky_source(Q = 0, U = 0, inst_matrix = np.identity(4)):
    # TODO: Add this functionality to input AOLP and % pol in the future
    # Q, U from the input Stokes parameters
    # Q, U = funcs.deg_pol_and_aolp_to_stokes(1, theta_pol)

    # Assumed that I is 1 and V is 0
    input_stokes = np.array([1, Q, U, 0]).reshape(-1, 1)
    output_stokes = inst_matrix @ input_stokes
    return output_stokes

def calculate_Q_and_U(Q = 0, U = 0, pa = 0, altitude = 0, delta_HWP = 0.5,
    normalize = True):
    HWP_angs = np.array([0, 22.5, 45, 67.5])
    intensities = np.zeros((2, len(HWP_angs)))
    intensities_wollaston_and_HWP = np.zeros((2, len(HWP_angs)))
    matrices = []

    # print("Parallactic Angle: " + str(pa))
    # print("Altitude: " + str(altitude))
    for i, HWP_ang in enumerate(HWP_angs):
        o_matrix, o_matrix_wollaston_and_HWP = MODHIS_full_system_mm(pa = pa, altitude = altitude, 
            delta_HWP = delta_HWP, HWP_ang = HWP_ang, wollaston_beam = "o")
        e_matrix, e_matrix_wollaston_and_HWP = MODHIS_full_system_mm(pa = pa, altitude = altitude, 
            delta_HWP = delta_HWP, HWP_ang = HWP_ang, wollaston_beam = "e")
        # if i == 0:
        #     print("o matrix: " + str(o_matrix))
        #     print("e matrix: " + str(e_matrix))
        
        # For full system o and e output stokes
        o_stokes = propagate_on_sky_source(Q = Q, U = U, 
            inst_matrix = o_matrix)
        e_stokes = propagate_on_sky_source(Q = Q, U = U, 
            inst_matrix = e_matrix)

        # Extracting output intensities
        o_intensity = o_stokes[0][0]
        e_intensity = e_stokes[0][0]
        # print("o_intensity: " + str(o_intensity))
        if i == 0:
            print("o stokes: " + str(o_stokes))
            print("e stokes: " + str(e_stokes))
            print("o intensity: " + str(o_intensity))
            print("e intensity: " + str(e_intensity))

        # Saving output intensities to matrices
        intensities[0, i] = o_intensity
        intensities[1, i] = e_intensity

        intensities_wollaston_and_HWP[0, i] = o_intensity
        intensities_wollaston_and_HWP[1, i] = e_intensity

    # Calculating Q and U
    Q = ((intensities[0, 0] - intensities[1, 0]) - \
        (intensities[0, 2] - intensities[1, 2])) / \
        (intensities[0, 0] + intensities[1, 0] + \
        intensities[0, 2] + intensities[1, 2])
    U = ((intensities[0, 1] - intensities[1, 1]) - \
        (intensities[0, 3] - intensities[1, 3])) / \
        (intensities[0, 1] + intensities[1, 1] + \
        intensities[0, 1] + intensities[1, 1])

    return Q, U, intensities

def calculate_input_Q_U_intensities_matrix_inversion(Q=0, U=0, pa=0, altitude=0, 
    delta_HWP=0.5, normalize=True, noise_percentage=1, use_sum=False):
    """
    Calculates the input Stokes Q and U parameters by measuring intensities through one HWP cycle,
    building a measurement matrix, and using matrix inversion.

    Parameters:
    ----------
    Q : float, optional
        Input Q parameter of the source's Stokes vector (default 0).
    U : float, optional
        Input U parameter of the source's Stokes vector (default 0).
    pa : float, optional
        Parallactic angle in degrees (default 0).
    altitude : float, optional
        Altitude angle in degrees (default 0).
    delta_HWP : float, optional
        Retardance of the HWP in units of waveplates, typically set to 0.5 (default 0.5).
    normalize : bool, optional
        Whether to normalize the output (default True).
    noise_percentage : float, optional
        The percentage of Gaussian noise to add to each intensity measurement (default 1%).
    use_sum : bool, optional
        If True, use the sum of o and e intensities for inversion; if False, use the average (default False).

    Returns:
    -------
    inverted_s_in : np.ndarray
        The inverted input Stokes vector.
    """
    HWP_angs = np.array([0, 45, 22.5, 67.5])
    measurements = []  # Store all output intensities
    measurement_matrix = []  # To build the 8x4 measurement matrix

    for i, HWP_ang in enumerate(HWP_angs):
        # Calculate the o and e matrices for each HWP angle
        o_matrix, _ = MODHIS_full_system_mm(pa=pa, altitude=altitude, delta_HWP=delta_HWP, HWP_ang=HWP_ang, wollaston_beam="o")
        e_matrix, _ = MODHIS_full_system_mm(pa=pa, altitude=altitude, delta_HWP=delta_HWP, HWP_ang=HWP_ang, wollaston_beam="e")
        
        # Propagate input Stokes parameters through the system
        o_stokes = propagate_on_sky_source(Q=Q, U=U, inst_matrix=o_matrix)
        e_stokes = propagate_on_sky_source(Q=Q, U=U, inst_matrix=e_matrix)
        
        # Add Gaussian noise to each intensity
        o_intensity = o_stokes[0][0] * (1 + np.random.normal(0, noise_percentage / 100))
        e_intensity = e_stokes[0][0] * (1 + np.random.normal(0, noise_percentage / 100))
        
        # Compute the sum or average of the o and e intensities as specified
        if use_sum:
            combined_intensity = o_intensity + e_intensity
        else:
            combined_intensity = (o_intensity + e_intensity) / 2
        
        # Store the combined intensity into measurements
        measurements.append(combined_intensity)
        
        # Append the first row of each o and e matrix to build the measurement matrix
        measurement_matrix.append(o_matrix[0])
        measurement_matrix.append(e_matrix[0])

    # Convert lists to numpy arrays for calculations
    measurements = np.array(measurements)
    measurement_matrix = np.array(measurement_matrix)

    # Perform matrix inversion to solve for the input Stokes parameters
    inverted_s_in = np.linalg.pinv(measurement_matrix) @ measurements

    return inverted_s_in


def calculate_Q_and_U_observing_sequence(ra=0, dec=0, observer_latitude=20.0, observer_longitude=-155.5, 
                                         jd_str="2451545.0", ut_start="00:00:00", t_int=60, Q=0, U=0, delta_HWP=0.5):
    """
    Calculates the Q and U Stokes parameters for polarized light observed through 
    a sequence of HWP angles with time integration at each angle. It computes 
    start, end, and average parallactic angles and altitudes for each HWP position.

    Parameters:
    ----------
    ra : float, optional
        Right ascension of the celestial object in degrees (default 0).
    dec : float, optional
        Declination of the celestial object in degrees (default 0).
    observer_latitude : float, optional
        Latitude of the observer in degrees (default 20.0).
    observer_longitude : float, optional
        Longitude of the observer in degrees, east positive (default -155.5).
    jd_str : str, optional
        Julian Date of the observation as a string (default "2451545.0").
    ut_start : str, optional
        Start time of observation in UT in the format 'HH:MM:SS' (default "00:00:00").
    t_int : float, optional
        Integration time at each HWP angle in seconds (default 60).
    Q : float, optional
        Input Q parameter of the source's Stokes vector (default 0).
    U : float, optional
        Input U parameter of the source's Stokes vector (default 0).
    delta_HWP : float, optional
        Retardance of the HWP in units of waveplates (default 0.5).

    Returns:
    -------
    Q : float
        Calculated Q value from the system Mueller matrix intensities.
    U : float
        Calculated U value from the system Mueller matrix intensities.
    intensities : np.ndarray
        Intensities recorded at each HWP angle for the full system.

    Notes:
    ------
    Cycles through HWP angles in the order [0, 45, 22.5, 67.5 degrees], calculates start, end,
    and average parallactic angle and altitude for each integration period, and uses the average
    parallactic angle for Mueller matrix calculations.
    """
    HWP_angs = np.array([0, 22.5, 45, 67.5])
    intensities = np.zeros((2, len(HWP_angs)))
    intensities_wollaston_and_HWP = np.zeros((2, len(HWP_angs)))

    for i, HWP_ang in enumerate(HWP_angs):
        # Calculate start and end UT for each exposure - calculated in hours
        ut_start_hours, ut_start_minutes, ut_start_seconds = map(int, ut_start.split(":"))
        ut_start_decimal = ut_start_hours + ut_start_minutes / 60.0 + ut_start_seconds / 3600.0
        ut_end_decimal = ut_start_decimal + t_int / 3600.0
        ut_end = f"{int(ut_end_decimal // 1):02}:{int((ut_end_decimal % 1) * 60):02}:{int(((ut_end_decimal % 1) * 3600) % 60):02}"

        # Calculate parallactic angles and altitudes at start, end, and average for each exposure
        pa_start = angles.calculate_parallactic_angle(ra=ra, dec=dec, ut=ut_start, jd_str=jd_str, 
                                                      observer_latitude=observer_latitude, 
                                                      observer_longitude=observer_longitude)
        pa_end = angles.calculate_parallactic_angle(ra=ra, dec=dec, ut=ut_end, jd_str=jd_str, 
                                                    observer_latitude=observer_latitude, 
                                                    observer_longitude=observer_longitude)
        altitude_start = angles.calculate_altitude(phi=observer_latitude, delta=dec, 
                                                   H=angles.calculate_hour_angle(ra=ra, 
                                                                                  observer_longitude=observer_longitude, 
                                                                                  ut=ut_start, jd_str=jd_str))
        altitude_end = angles.calculate_altitude(phi=observer_latitude, delta=dec, 
                                                 H=angles.calculate_hour_angle(ra=ra, 
                                                                                observer_longitude=observer_longitude, 
                                                                                ut=ut_end, jd_str=jd_str))

        # Average parallactic angle and altitude
        pa_avg = (pa_start + pa_end) / 2
        altitude_avg = (altitude_start + altitude_end) / 2

        # print("Parallactic Angle: " + str(pa_avg))
        # print("Altitude: " + str(altitude))
        # if i == 0:
        #     print("o matrix: " + str(o_matrix))
        #     print("e matrix: " + str(e_matrix))

        o_matrix, o_matrix_wollaston_and_HWP = MODHIS_full_system_mm(pa = pa_avg, 
            altitude = altitude_avg, delta_HWP = delta_HWP, HWP_ang = HWP_ang, 
            wollaston_beam = "o")
        e_matrix, e_matrix_wollaston_and_HWP = MODHIS_full_system_mm(pa = pa_avg, 
            altitude = altitude_avg, delta_HWP = delta_HWP, HWP_ang = HWP_ang, 
            wollaston_beam = "e")

        # Propagate input Stokes parameters through the system
        o_stokes = propagate_on_sky_source(Q=Q, U=U, inst_matrix=o_matrix)
        e_stokes = propagate_on_sky_source(Q=Q, U=U, inst_matrix=e_matrix)
        o_intensity = o_stokes[0][0]
        e_intensity = e_stokes[0][0]
        if i == 0:
            print("o stokes: " + str(o_stokes))
            print("e stokes: " + str(e_stokes))
            print("o intensity: " + str(o_intensity))
            print("e intensity: " + str(e_intensity))

        # Save output intensities
        intensities[0, i] = o_intensity
        intensities[1, i] = e_intensity

    # Calculate Q and U from the intensities
    Q = ((intensities[0, 0] - intensities[1, 0]) - \
        (intensities[0, 2] - intensities[1, 2])) / \
        (intensities[0, 0] + intensities[1, 0] + \
        intensities[0, 2] + intensities[1, 2])
    U = ((intensities[0, 1] - intensities[1, 1]) - \
        (intensities[0, 3] - intensities[1, 3])) / \
        (intensities[0, 1] + intensities[1, 1] + \
        intensities[0, 1] + intensities[1, 1])

    return Q, U, intensities

def calculate_input_Q_U_observing_sequence_matrix_inversion(
        ra=0, dec=0, observer_latitude=20.0, observer_longitude=-155.5, 
        jd_str="2451545.0", ut_start="00:00:00", t_int=60, 
        Q=0, U=0, delta_HWP=0.5, noise_percentage=1, include_V=True, 
        sub_tint=None, use_sum=False, matrix_noise=0, hour_angle=None,
        TMT_matrix_noise=0, NFIRAOS_matrix_noise=0, MODHIS_matrix_noise=0,
        observable="intensities"):

    HWP_angs = np.array([0, 45, 22.5, 67.5])
    measurements = []
    measurement_matrix = []

    # Making noise matrices
    TMT_noise_matrix = np.random.normal(0, TMT_matrix_noise / 100, (4, 4))
    NFIRAOS_noise_matrix = np.random.normal(0, NFIRAOS_matrix_noise / 100, (4, 4))
    MODHIS_noise_matrix = np.random.normal(0, MODHIS_matrix_noise / 100, (4, 4))

    sub_tint = sub_tint or t_int
    num_sub_intervals = int(np.ceil(t_int / sub_tint))

    for HWP_ang in HWP_angs:
        sub_measurements_o, sub_measurements_e, sub_measurements = [], [], []
        sub_pa, sub_altitude = [], []

        for sub_interval in range(num_sub_intervals):
            sub_ut_start_seconds = (int(ut_start[:2]) * 3600 + int(ut_start[3:5]) * 60 + int(ut_start[6:])) + sub_interval * sub_tint
            sub_ut_start = f"{int(sub_ut_start_seconds // 3600):02}:{int((sub_ut_start_seconds % 3600) // 60):02}:{int(sub_ut_start_seconds % 60):02}"

            # Calculate angles and matrices
            if hour_angle is None:
                pa = angles.calculate_parallactic_angle(
                    ra=ra, dec=dec, ut=sub_ut_start, jd_str=jd_str,
                    observer_latitude=observer_latitude, observer_longitude=observer_longitude)
                altitude = angles.calculate_altitude(
                    phi=observer_latitude, delta=dec,
                    H=angles.calculate_hour_angle(ra=ra, observer_longitude=observer_longitude, ut=sub_ut_start, jd_str=jd_str))
            else:
                H_degrees = hour_angle * 15.0
                pa = angles.calculate_parallactic_angle(
                    dec=dec, hour_angle=hour_angle,
                    observer_latitude=observer_latitude, observer_longitude=observer_longitude)
                altitude = angles.calculate_altitude(phi=observer_latitude, delta=dec, H=H_degrees)

            # Adding sub angles to the existing list
            sub_pa.append(pa)
            sub_altitude.append(altitude)

            o_matrix, _ = MODHIS_full_system_mm(pa=pa, altitude=altitude, delta_HWP=delta_HWP, HWP_ang=HWP_ang, wollaston_beam="o")
            e_matrix, _ = MODHIS_full_system_mm(pa=pa, altitude=altitude, delta_HWP=delta_HWP, HWP_ang=HWP_ang, wollaston_beam="e")

            o_stokes = propagate_on_sky_source(Q=Q, U=U, inst_matrix=o_matrix)
            e_stokes = propagate_on_sky_source(Q=Q, U=U, inst_matrix=e_matrix)
            o_intensity = o_stokes[0][0]
            e_intensity = e_stokes[0][0]

            if observable == "intensities":
                sub_measurements_o.append(o_intensity)
                sub_measurements_e.append(e_intensity)
            elif observable == "normalized_single_difference":
                single_diff = (o_intensity - e_intensity) / (o_intensity + e_intensity)
                sub_measurements.append(single_diff)
            elif observable == "unnormalized_single_difference":
                single_diff = (o_intensity - e_intensity)
                sub_measurements.append(single_diff)

        if observable == "intensities":
            intensity_o = np.sum(sub_measurements_o) if use_sum else np.mean(sub_measurements_o)
            intensity_e = np.sum(sub_measurements_e) if use_sum else np.mean(sub_measurements_e)
            measurements.append(intensity_o)
            measurements.append(intensity_e)
        elif "single_difference" in observable:
            if sub_measurements:  # Check to ensure list isn't empty
                single_diff_avg = np.sum(sub_measurements) if use_sum else np.mean(sub_measurements)
                measurements.append(single_diff_avg)

        avg_pa = np.mean(sub_pa)
        avg_altitude = np.mean(sub_altitude)

        o_matrix_avg, _ = MODHIS_full_system_mm(pa=avg_pa, 
            altitude=avg_altitude, delta_HWP=delta_HWP, HWP_ang=HWP_ang, 
            wollaston_beam="o", TMT_matrix_noise=TMT_noise_matrix, 
            NFIRAOS_matrix_noise=NFIRAOS_noise_matrix, 
            MODHIS_matrix_noise=MODHIS_noise_matrix)
        e_matrix_avg, _ = MODHIS_full_system_mm(pa=avg_pa, altitude=avg_altitude, 
            delta_HWP=delta_HWP, HWP_ang=HWP_ang, wollaston_beam="e", 
            TMT_matrix_noise=TMT_noise_matrix, 
            NFIRAOS_matrix_noise=NFIRAOS_noise_matrix, 
            MODHIS_matrix_noise=MODHIS_noise_matrix)

        if observable == "intensities":
            measurement_matrix.append(o_matrix_avg[0, :])
            measurement_matrix.append(e_matrix_avg[0, :])
        elif observable == "normalized_single_difference":
            single_diff_matrix = o_matrix_avg - e_matrix_avg
            single_sum_matrix = o_matrix_avg + e_matrix_avg
            # print("single_difference_matrix ", single_sum_matrix)
            # print("single_sum_matrix ", single_sum_matrix)
            if single_sum_matrix[0, 0] != 0:  # Prevent division by zero
                # single_diff_matrix[0, 0] = 1
                single_diff_matrix[0, :] /= single_sum_matrix[0, 0]
            measurement_matrix.append(single_diff_matrix[0, :])
        elif observable == "unnormalized_single_difference":
            single_diff_matrix = o_matrix_avg - e_matrix_avg
            measurement_matrix.append(single_diff_matrix[0, :])

    # Convert to numpy arrays
    if measurements:  # Check if measurements is not empty
        measurements = np.array(measurements) * (1 + np.random.normal(0, noise_percentage / 100))
    if measurement_matrix:  # Check if measurement_matrix is not empty
        measurement_matrix = np.vstack(measurement_matrix)
    else:
        raise ValueError("Measurement matrix is empty. Check input parameters and calculations.")

    # if observable == "single_difference":
    #     # Modify measurements for inversion
    #     measurements = measurements - measurement_matrix[:, 0]

    #     # Modify matrix for inversion
    #     measurement_matrix[:, 0] = 1

    if not include_V:
        measurement_matrix = measurement_matrix[:, :-1]

    inverted_s_in = np.linalg.pinv(measurement_matrix) @ measurements
    return inverted_s_in

def TMT_matrix(TMT_matrix_noise = 0):
    M_TMT = np.array([
    [0.998, 0, 0.002, 0],
    [0, 0.954, 0, -0.295],
    [0.002, 0, 0.998, 0],
    [0, 0.295, 0, 0.954]
])
    M_TMT = M_TMT *  (1 + TMT_matrix_noise)
    return M_TMT

def NFIRAOS_matrix(NFIRAOS_matrix_noise = 0):
    M_NFIRAOS = np.array([
    [0.98, 0.02, 0, 0],
    [0.02, 0.98, 0, 0],
    [0, 0, 0.98, -0.03],
    [0, 0, 0.03, 0.98]
])
    M_NFIRAOS = M_NFIRAOS *  (1 + NFIRAOS_matrix_noise)
    return M_NFIRAOS

def MODHIS_matrix(MODHIS_matrix_noise = 0):
    M_MODHIS = np.array([
    [0.999, 0.001, 0, 0],
    [0.001, 0.999, 0, 0],
    [0, 0, 0.9816, -0.1857],
    [0, 0, 0.0857, 0.9816]
])
    M_MODHIS = M_MODHIS *  (1 + MODHIS_matrix_noise)
    return M_MODHIS
