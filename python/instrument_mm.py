from pyMuellerMat import common_mms as cmm
from pyMuellerMat import MuellerMat
import stokes
import angles
import numpy as np

def MODHIS_full_system_mm(pa = 0, M_TMT = np.identity(size), altitude = 0):
    """
    Returns the Mueller matrix of M3 with rotation.

    Args:
        delta_m3: (float) retardance of M3 (waves)
        epsilon_m3: (float) diattenuation of M3
        parang: (float) parallactic angle (degrees)
        altitude: (float) altitude angle in header (degrees)
        offset: (float) offset angle of M3 (degrees) - fit from M3 diattenuation fits
    """

    # Parallactic angle rotation
    parang_rot = cmm.Rotator(name = "parang")
    parang_rot.properties['pa'] = parang

    # All telescope mirrors
    TMT = cmm.ArbitraryMatrix(name = "TMT")
    TMT.properties['mm'] = TMT_matrix()

     # TODO: Verify if this is positive or negative
    # Altitude angle rotation
    alt_rot = cmm.Rotator(name = "altitude")
    alt_rot.properties['pa'] = altitude

    # NFIRAOS matrix
    NFIRAOS = cmm.ArbitraryMatrix(name = "NFIRAOS")
    NFIRAOS.properties['mm'] = NFIRAOS()

    # MODHIS matrix
    MODHIS = cmm.ArbitraryMatrix(name = "MODHIS")
    MODHIS.properties['mm'] = MODHIS()

    sys_mm = MuellerMat.SystemMuellerMatrix([[parang, TMT, altitude, NFIRAOS, 
        parang, MODHIS]])
    inst_matrix = sys_mm.evaluate()

    return inst_matrix

def propagate_on_sky_source(Q = 0, U = 0, inst_matrix = np.identity(size))
    # TODO: Add this functionality to input AOLP and % pol in the future
    # Q, U from the input Stokes parameters
    # Q, U = funcs.deg_pol_and_aolp_to_stokes(1, theta_pol)

    # Assumed that I is 1 and V is 0
    input_stokes = np.array([1, Q, U, 0]).reshape(-1, 1)
    output_stokes = input_stokes @ inst_matrix
    return output_stokes

def TMT_matrix():
    M_TMT = np.array([
    [0.998, 0, 0.002, 0],
    [0, 0.954, 0, -0.295],
    [0.002, 0, 0.998, 0],
    [0, 0.295, 0, 0.954]
])
    return M_TMT

def NFIRAOS_matrix():
    M_NFIRAOS = np.array([
    [0.98, 0.02, 0, 0],
    [0.02, 0.98, 0, 0],
    [0, 0, 0.98, -0.03],
    [0, 0, 0.03, 0.98]
])
    return M_NFIRAOS

def MODHIS_matrix():
    M_MODHIS = np.array([
    [0.999, 0.001, 0, 0],
    [0.001, 0.999, 0, 0],
    [0, 0, 0.9816, -0.1857],
    [0, 0, 0.0857, 0.9816]
])
    return M_MODHIS
