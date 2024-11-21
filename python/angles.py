import math

def calculate_altitude(phi=20.0, delta=0.0, H=0.0):
    """
    Calculate the altitude (E) of a celestial object given the observer's latitude (phi),
    the object's declination (delta), and the hour angle (H).

    Parameters:
    ----------
    phi : float, optional
        Latitude of the observer in degrees (default is 20.0, approx. Mauna Kea).
    delta : float, optional
        Declination of the target object in degrees (default is 0.0, equatorial object).
    H : float, optional
        Hour angle of the target object in degrees (default is 0.0, transit time).

    Returns:
    -------
    E : float
        Altitude of the target object in degrees.
    """
    phi_rad = math.radians(phi)
    delta_rad = math.radians(delta)
    H_rad = math.radians(H)
    
    sin_E = math.sin(phi_rad) * math.sin(delta_rad) + math.cos(phi_rad) * math.cos(delta_rad) * math.cos(H_rad)
    sin_E = min(1.0, max(-1.0, sin_E))  # Clamp value to avoid errors
    E_deg = math.degrees(math.asin(sin_E))
    
    return E_deg

def calculate_hour_angle(ra=0.0, observer_longitude=-155.5, ut="00:00:00", jd_str="2451545.0"):
    ra_hours = ra / 15.0
    ut_hours, ut_minutes, ut_seconds = map(int, ut.split(":"))
    ut_decimal = ut_hours + ut_minutes / 60 + ut_seconds / 3600
    jd = float(jd_str)
    
    # Calculate Julian centuries since J2000.0
    T = (jd - 2451545.0) / 36525.0
    
    # Calculate GMST at 0h UT using more precise formula
    gmst_0h = 100.46061837 + 36000.770536 * T + 0.000387933 * T**2 - T**3 / 38710000.0
    gmst_0h = gmst_0h % 360  # Ensure it is within 0-360 degrees
    
    # Calculate GMST at UT
    gmst = gmst_0h + 360.98564736629 * (ut_decimal / 24.0)
    gmst = gmst % 360  # Ensure it is within 0-360 degrees
    gmst_hours = gmst / 15.0  # Convert to hours
    
    # Calculate Local Sidereal Time (LST)
    lst = (gmst_hours + observer_longitude / 15.0) % 24
    
    # Calculate Hour Angle (H)
    H = (lst - ra_hours) % 24
    H_degrees = H * 15.0
    
    return H_degrees

def calculate_parallactic_angle(ra=0.0, dec=0.0, ut="00:00:00", jd_str="2451545.0", observer_latitude=20.0, observer_longitude=-155.5, hour_angle = None):
    """
    Calculate the parallactic angle of a celestial object given its right ascension (RA),
    declination (Dec), observer's latitude, longitude, UT, and Julian Date (JD).

    Parameters:
    ----------
    ra : float, optional
        Right ascension of the celestial object in degrees (default is 0.0).
    dec : float, optional
        Declination of the celestial object in degrees (default is 0.0).
    ut : str, optional
        Time in UT as 'HH:MM:SS' (default is "00:00:00").
    jd_str : str, optional
        Julian Date as a string (default is "2451545.0", J2000).
    observer_latitude : float, optional
        Latitude of the observer in degrees (default is 20.0, approx. Mauna Kea).
    observer_longitude : float, optional
        Longitude of the observer in degrees, east positive (default is -155.5, Mauna Kea).

    Returns:
    -------
    q : float
        Parallactic angle of the celestial object in degrees.
    """
    if hour_angle == None:
        H_degrees = calculate_hour_angle(ra=ra, observer_longitude=observer_longitude, ut=ut, jd_str=jd_str)
    else:
        H_degrees = hour_angle * 15.0
    H_rad = math.radians(H_degrees)
    lat_rad = math.radians(observer_latitude)
    dec_rad = math.radians(dec)
    
    tan_q = math.sin(H_rad) / (math.tan(lat_rad) * math.cos(dec_rad) - math.sin(dec_rad) * math.cos(H_rad))
    q_deg = math.degrees(math.atan(tan_q))
    
    return q_deg
