def stokes_to_deg_pol_and_aolp(Q, U):
    pol_percent = np.sqrt(Q ** 2 + U ** 2) * 100  # Convert to percentage
    aolp = 0.5 * np.arctan2(U, Q) * (180/np.pi)  # Convert to degrees
    return pol_percent, aolp

def deg_pol_and_aolp_to_stokes(pol_percent, aolp):
    # Convert percentage polarization to a fraction
    pol_fraction = pol_percent / 100.0
    
    # Convert aolp from degrees to radians
    aolp_rad = np.deg2rad(aolp * 2)  # Factor of 2 due to the 0.5 factor in arctan2

    # Calculate Q and U
    Q = pol_fraction * np.cos(aolp_rad)
    U = pol_fraction * np.sin(aolp_rad)

    return Q, U