import math

def calculate_altitude(phi, delta, H):
    """
    Calculate the altitude (E) of a celestial object given the observer's latitude (phi), 
    the object's declination (delta), and the hour angle (H).
    
    Parameters:
    ----------
    phi : float
        The latitude of the observer, in degrees.
    delta : float
        The declination of the target object, in degrees.
    H : float
        The hour angle of the target object, in degrees.
    
    Returns:
    -------
    E : float
        The altitude of the target object, in degrees.
    """
    
    # Convert degrees to radians for trigonometric calculations
    phi_rad = math.radians(phi)
    delta_rad = math.radians(delta)
    H_rad = math.radians(H)
    
    # Calculate sin(E) using the given formula
    sin_E = math.sin(phi_rad) * math.sin(delta_rad) + math.cos(phi_rad) * math.cos(delta_rad) * math.cos(H_rad)
    
    # Ensure the value is in the range [-1, 1] to avoid numerical errors with asin
    sin_E = min(1.0, max(-1.0, sin_E))
    
    # Calculate the altitude (E) in radians and convert to degrees
    E_rad = math.asin(sin_E)
    E_deg = math.degrees(E_rad)
    
    return E_deg

def calculate_hour_angle(ra, observer_longitude, ut, jd_str):
    """
    Calculate the hour angle (H) of a celestial object given its right ascension (RA), 
    the observer's longitude, the time of observation in UT, and the Julian Date (JD).

    Parameters:
    ----------
    ra : float
        The right ascension of the celestial object, in degrees.
    observer_longitude : float
        The longitude of the observer, in degrees (positive for east, negative for west).
    ut : str
        The time of observation in UT as a string in the format 'HH:MM:SS'.
    jd_str : str
        The Julian Date of the observation as a string.

    Returns:
    -------
    H : float
        The hour angle of the celestial object, in degrees.
    """

    # Convert RA from degrees to hours
    ra_hours = ra / 15.0

    # Parse UT string to get hours, minutes, and seconds
    ut_hours, ut_minutes, ut_seconds = [int(part) for part in ut.split(":")]
    ut_decimal = ut_hours + ut_minutes / 60 + ut_seconds / 3600

    # Convert Julian Date from string to float
    jd = float(jd_str)

    # Calculate the number of days since J2000.0 (JD for J2000.0 is 2451545.0)
    jd_j2000 = jd - 2451545.0

    # Calculate the Greenwich Sidereal Time (GST) in degrees at 0h UT
    gst_0h = 100.46061837 + 36000.770053608 * (jd_j2000 / 36525.0)
    gst_0h = gst_0h % 360  # Normalize to 0-360 degrees

    # Update GST to the observation time in UT
    gst = (gst_0h + 360.98564724 * (ut_decimal / 24.0)) % 360

    # Convert GST from degrees to hours
    gst_hours = gst / 15.0

    # Calculate Local Sidereal Time (LST)
    lst = (gst_hours + observer_longitude / 15.0) % 24

    # Calculate the hour angle (in hours) and convert to degrees
    H = (lst - ra_hours) % 24
    H_degrees = H * 15.0

    return H_degrees

def calculate_parallactic_angle(ra, dec, ut, jd_str, observer_latitude, observer_longitude):
    """
    Calculate the parallactic angle of a celestial object given its right ascension (RA), 
    declination (Dec), the observer's latitude, longitude, and the time of observation in UT 
    and Julian Date (JD).

    Parameters:
    ----------
    ra : float
        The right ascension of the celestial object, in degrees.
    dec : float
        The declination of the celestial object, in degrees.
    ut : str
        The time of observation in UT as a string in the format 'HH:MM:SS'.
    jd_str : str
        The Julian Date of the observation as a string.
    observer_latitude : float
        The latitude of the observer, in degrees.
    observer_longitude : float
        The longitude of the observer, in degrees (positive for east, negative for west).

    Returns:
    -------
    q : float
        The parallactic angle of the celestial object, in degrees.
    """
    
    # Convert RA and Dec from degrees to radians
    ra_rad = math.radians(ra)
    dec_rad = math.radians(dec)
    
    # Calculate the Hour Angle (H) using the previous hour angle function logic
    H_degrees = calculate_hour_angle(ra, observer_longitude, ut, jd_str)
    H_rad = math.radians(H_degrees)
    
    # Convert observer's latitude to radians
    lat_rad = math.radians(observer_latitude)
    
    # Calculate the parallactic angle using the formula
    sin_H = math.sin(H_rad)
    cos_H = math.cos(H_rad)
    tan_phi = math.tan(lat_rad)
    cos_delta = math.cos(dec_rad)
    sin_delta = math.sin(dec_rad)
    
    # Parallactic angle formula
    tan_q = sin_H / (tan_phi * cos_delta - sin_delta * cos_H)
    q_rad = math.atan(tan_q)
    
    # Convert the parallactic angle to degrees
    q_deg = math.degrees(q_rad)
    
    return q_deg

# You can use the calculate_hour_angle function provided earlier in the same script

