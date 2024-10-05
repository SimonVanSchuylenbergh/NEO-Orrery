import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import functools
from pandas import to_datetime
from itertools import combinations, chain
from sys import float_info
from math import isinf

# WMPG was here

#Equations for parabolic and hyperbolic orbits taken from
#https://www.bogan.ca/orbits/kepler/orbteqtn.html
#Website source reference is:
#Reference: "Fundamentals of Astrodynamics"
#by R.R.Bate, D.D.Mueller, and J.E. White
#Dover Publications (1971)

np.set_printoptions(suppress=True)
neg_zero_eps = -1e-14 #negative values smaller than this are set to 0
zero_eps = 1e-14 #positive values smaller than this are considered to be zero

#Places the arrays of observational data into a dictionary and sorts each object's observations by Julian Date
#Output: asteroids (dictionary of observational data with each entry being its own object)
def initialize_obs_data():
    #5 observations each --> (Julian Date, Right Ascension, Declination, Magnitude)
    X0247 = np.array([
    [2452465.5, 93.9350137072466680, 38.1582374014492629, 12.026],
    [2452470.5, 95.6958342344323825, 38.0852446940591847, 12.043],
    [2452480.5, 99.0739551229480071, 37.9339616031920031, 12.063],
    [2452468.5, 94.9970501725531307, 38.1147826741125755, 12.037],
    [2452474.5, 97.0706479192634788, 38.0251849185151869, 12.053]], dtype=np.float64)

    X0487 = np.array([
    [2452465.5, 259.7657415651311794, -12.6153217571698129, 18.047],
    [2452470.5, 259.2650934782980130, -12.6523752347416423, 18.059],
    [2452480.5, 258.3754093732495676, -12.7520580832093593, 18.085],
    [2452468.5, 259.4612764378233010, -12.6365038571287780, 18.054],
    [2452474.5, 258.8902468513585973, -12.6882436206546370, 18.069]], dtype=np.float64)

    X0506 = np.array([
    [2452465.5, 1.3498908499056173, -17.0544601237825759, 34.232],
    [2452470.5, 1.3471518603560630, -17.0588726579291183, 34.231],
    [2452480.5, 1.3392532114351825, -17.0683461604096891, 34.229],
    [2452468.5, 1.3483466887310642, -17.0570786799518146, 34.231],
    [2452474.5, 1.3443720077316741, -17.0625675299641095, 34.230]], dtype=np.float64)

    X0613 = np.array([ #this was the object specifically assigned to me
    [2452465.5, 167.0760309907474266, 1.1312156781805742, 23.377],
    [2452470.5, 167.1707119127673593, 1.1081564268121049, 23.370],
    [2452480.5, 167.3941338167601032, 1.0452617836226097, 23.353],
    [2452468.5, 167.1314041099444125, 1.1180749225401012, 23.373],
    [2452474.5, 167.2548550540942358, 1.0856089445055799, 23.364]], dtype=np.float64)

    X0718 = np.array([
    [2452465.5, 15.6012629718496640, -9.4957892151581405, 25.010],
    [2452470.5, 15.6120191227987597, -9.5276575325865682, 25.006],
    [2452480.5, 15.6078224720019083, -9.6009907311125691, 24.996],
    [2452468.5, 15.6087463680267593, -9.5145048425643370, 25.008],
    [2452474.5, 15.6144447236321664, -9.5555221324367370, 25.003]], dtype=np.float64)

    X0769 = np.array([
    [2452465.5, 357.6833615090893659, 0.3618335233054201, 11.951],
    [2452470.5, 357.8282418937238845, 0.4280983288767462, 11.948],
    [2452480.5, 357.7820777378383355, 0.4111137264269742, 11.928],
    [2452468.5, 357.7834543534117984, 0.4074229619432249, 11.950],
    [2452474.5, 357.8642453789895512, 0.4456369974434286, 11.942]], dtype=np.float64)

    asteroids = dict([('X0247', X0247), ('X0487', X0487), ('X0506', X0506), ('X0613', X0613), ('X0718', X0718), ('X0769', X0769)])

    for data in asteroids.values(): data.view('f8,f8,f8,f8').sort(order=['f0'], axis=0) #sort the arrays in place by their Julian Date

    return asteroids

#Creates a dictionary of the orbital elements for each of the planets (plus Pluto), then modifies it to be in a more useful form
#Output: planets (the dictionary of orbital data, plus the name of the object, and its plotted orbit color)
#Output: epoch (the Julian Date associated with the provided orbital elements)
def initialize_planet_data():
    #dictionary data --> key: name, values: ploted orbit color, orbital element data (numpy array)
    #numpy data --> semi-major axis (AU), eccentricity, inclination (deg), longitude of the ascending node (deg), longitude of perihelion (deg), mean longitude (deg)

    epoch = 2451545 #Julian Date for Jan 1, 2000 at 12h (J2000 epoch)

    planets = dict([
    ("Mercury", ["dimgrey", np.array([0.38709893, 0.20563069, 7.00487, 48.33167, 77.45645, 252.25084])]),
    ("Venus", ["chocolate", np.array([0.72333199, 0.00677323, 3.39471, 76.68069, 131.53298, 181.97973])]),
    ("Earth", ["mediumblue", np.array([1.00000011, 0.01671022, 0.00005, -11.26064, 102.94719, 100.46435])]),
    ("Mars", ["firebrick", np.array([1.52366231, 0.09341233, 1.85061, 49.57854, 336.04084, 355.45332])]),
    ("Jupiter", ["darkorange", np.array([5.20336301, 0.04839266, 1.30530, 100.55615, 14.75385, 34.40438])]),
    ("Saturn", ["peachpuff", np.array([9.53707032, 0.05415060, 2.48446, 113.71504, 92.43194, 49.94432])]),
    ("Uranus", ["lightblue", np.array([19.19126393, 0.04716771, 0.76986, 74.22988, 170.96424, 313.23218])]),
    ("Neptune", ["royalblue", np.array([30.06896348, 0.00858587, 1.76917, 131.72169, 44.97135, 304.88003])]),
    ("Pluto", ["tan", np.array([39.48168677, 0.24880766, 17.14175, 110.30347, 224.06676, 238.92881])])])

    for _, data in planets.values():
        data[2:] *= np.pi/180 #convert all degrees to radians
        data[-1] -= data[-2] #change the mean longitude to be the mean anomaly (L = b_omega + s_omega + M --> M)
        data[-2] -= data[-3] #change the longitude of perihelion to be the argument of perihelion (s_omega_hat = b_omega + s_omega --> s_omega)

    #numpy data is now --> semi-major axis (AU), eccentricity, inclination (rad), longitude of the ascending node (rad), argument of perihelion (rad), mean anomaly (rad)
    return planets, epoch

#Computes the Right Ascension and Declination of an object, given its Cartsian coordinates (from assignment 2)
#Input: (x, y, z) Cartsian coordinates in the geocentric equatorial / ecliptic coordinate system
#Output: Right Ascension (degrees) / ecliptic latitude (degrees)
#Output: Declination (degrees) / ecliptic longitude (degrees)
def ra_and_dec(x, y, z):
    r = np.sqrt(x*x + y*y + z*z)
    dec = np.arcsin(z / r) #can represent ecliptic latitude 
    ra = np.arctan2(y, x) #can represent ecliptic longitude

    dec *= 180 / np.pi #convert to degrees
    ra *= 180 / np.pi #convert to degrees

    if hasattr(dec, '__len__'):
        dec[dec > 90] -= 360 #ensure Dec is within -90 to 90 degrees
        ra[ra < 0] += 360 #RA needs to be converted to all positive bounds
    else: #not an array
        if dec > 90: dec -= 360
        if ra < 0: ra += 360

    return (ra, dec)

#Computes the the eccentric anomaly, given the true anomaly along with the eccentricity (modified form assignment 5)
#All input angles are in radians
#Input: f (the true anomaly --> units: radians)
#Input: e (eccentricity)
#Input: force_positive (boolean toggle to force the output to be positive --> 0 to 2*pi)
#Output: E (the eccentric anomaly --> units: radians)
def compute_eccentric_anomaly(f, e, force_positive=False):
    if e < 1: #elliptical orbit
        E = 2*np.arctan(np.sqrt((1-e) / (1+e)) * np.tan(f/2))
        
        if force_positive: #ensure positive valued angles
            if hasattr(E, '__len__'): E[E < 0] += 2*np.pi #array
            else: E += 2*np.pi * (E < 0) #single value
    elif e > 1: E = 2*np.arctanh(np.sqrt((e-1) / (1+e)) * np.tan(f/2)) #hyperbolic orbit
    else: E = np.tan(f/2) #parabolic orbit
    
    return E

#Computes the the true anomaly, given the eccentric anomaly along with the eccentricity (modified form assignment 5)
#All input angles are in radians
#Input: E (the eccentric anomaly --> units: radians)
#Input: e (eccentricity)
#Input: force_positive (boolean toggle to force the output to be positive --> 0 to 2*pi)
#Output: f (the true anomaly --> units: radians)
def compute_true_anomaly(E, e, force_positive=False):
    if e < 1: f = 2*np.arctan(np.sqrt((1+e) / (1-e)) * np.tan(E/2)) #elliptical orbit
    elif e > 1: f = 2*np.arctan(np.sqrt((e+1) / (e-1))*np.tanh(E/2)) #hyperbolic orbit
    else: f = 2*np.arctan(E) #parabolic orbit

    if force_positive: #ensure positive valued angles
        if hasattr(f, '__len__'): f[f < 0] += 2*np.pi #array
        else: f += 2*np.pi * (f < 0) #single value
    
    return f

#Computes the standard orbital elements of an object given its heliocentric ecliptic position and velocity vectors (from assignment 4)
#along with mu = G(m1 + m2). All inputs and outputs are in astronomically useful units:
#length is in AU, mass is in solar masses, and time is in days. G is therefore Gauss’ gravitational constant squared
#Input: mu (Gauss’ gravitational constant squared times the sum of the masses of the orbiting bodies)
#Input: x, y, z (the 3 elements of the position vector, in units of AU)
#Input: vx, vy, vz (the 3 elements of the velocity vector, in units of AU/day)
#Output: a (semi-major axis, in AU)
#Output: e (eccentricity, unitless)
#Output: i (inclination, degrees)
#Output: b_omega (longitude of the ascending node, degrees)
#Output: s_omega (argument of periapsis, degrees)
#Output: f (true anomaly, degrees)
def compute_orbital_elements(mu, x, y, z, vx, vy, vz):
    r = np.sqrt(x*x + y*y + z*z) #distance
    v2 = vx*vx + vy*vy + vz*vz #velocity squared
    rdv = x*vx + y*vy + z*vz #r dot v
    
    #angular momentum per unit mass
    hx = y*vz - z*vy
    hy = z*vx - x*vz
    hz = x*vy - y*vx

    #intermediate momentum related values
    hxy2 = hx*hx + hy*hy
    hxy = np.sqrt(hxy2)
    h2 = hxy2 + hz*hz
    h = np.sqrt(h2)

    if np.abs(r) < zero_eps: #Inside the Sun
        print("It is just me or is it hot in here?") #what is this an orbit for ants?
        return (None,)*6

    v2_crit = 2*mu/r #the velocity squared at which the orbit becomes parabolic (above which it is hyperbolic)

    if np.abs(h) < zero_eps: #not an orbit, but rather a straight line
        if v2 - v2_crit > neg_zero_eps and rdv > 0: print("Orbital collapse, escaping system")
        else: print("Orbital collapse, collision inevitable")
        return (None,)*6

    if np.abs(v2 - v2_crit) < zero_eps: #parabolic orbit:
        print("Warning: This object is on a parabolic orbit, and will return a value of inf for the semi-major axis")
        a = float('inf')
        e = 1 #a parabolic orbit has an eccentricity of exactly 1
    else:
        a = 1 / (2/r - v2/mu) #semi-major axis
        e2 = 1 - h2/(mu*a) #eccentricity squared
        if e2 < 0: e2 = 0 #ensures e2 is greater than zero, can be less than zero due to floating point rounding errors (e2 *= e2 > 0 sets it to negative 0 for some reason)
        e = np.sqrt(e2) #eccentricity (can be greater than 1 --> hyperbolic)
    
    if v2 - v2_crit > zero_eps:
        print("Warning: This object is on a hyperbolic orbit, and will return a negative value for the semi-major axis")

    i = np.rad2deg(np.arctan2(hxy, hz)) #inclination (range 0 to 180 degrees)

    if i == 0 or i == 180:
        b_omega = 0 #set the longitude of the ascending node to zero when the inclination is zero
        u = np.arccos(x/r)
        u *= np.sign(y) #u is now in the same range as the output of arctan2 (-pi to pi) {x | (y & (1 << 31)) fast c++ version}
        if i == 180: u *= -1 #retrograde flip
    else:
        b_omega = np.arctan2(hx, -hy) #longitude of the ascending node
        u = np.arctan2(z*h/hxy, x*np.cos(b_omega) + y*np.sin(b_omega)) #argument of latitude (argument of periapsis + true anomaly)

        b_omega = np.rad2deg(b_omega) #convert from radians to degrees
        if neg_zero_eps < b_omega < 0: b_omega = 0 #fix floating point issues
        b_omega += 360 * (b_omega < 0) #set within the range of 0 to 360 degrees

    if e == 0:
        s_omega = 0 #in a circular orbit there is no periapsis, therefore the argument of periapsis is set to zero
        f = u #the true anomaly is therefore just the argument of latitude
    else:
        f = np.arctan2(rdv*h, h2 - mu*r) #true anomaly
        s_omega = u - f #argument of periapsis

        s_omega = np.rad2deg(s_omega) #convert from radians to degrees
        if neg_zero_eps < s_omega < 0: s_omega = 0 #fix floating point issues
        s_omega += 360 * (s_omega < 0) #set within the range of 0 to 360 degrees
        s_omega -= 360 * (s_omega == 360) #extremely rare but possible for this to occur due to u - f (set 360 to 0)

    f = np.rad2deg(f) #convert from radians to degrees
    if neg_zero_eps < f < 0: f = 0 #fix floating point issues
    f += 360 * (f < 0) #set within the range of 0 to 360 degrees

    return (a, e, i, b_omega, s_omega, f)

#Computes the postion and (optionally) velocity of an orbiting object around a central mass (e.g Earth around the Sun) in the
#orbital plane, given the object's orbital elements a, e, and E, along with the standard gravitational parameter of the central mass (modified from assignment 5)
#Input: mu (standard gravitational parameter of the central mass, or object-central mass system)
#Input: a (semi-major axis, units: AU)
#Input: e (eccentricity)
#Input: E (eccentric anomaly, units: radians)
#Input: f (true anomaly, units: radians)
#Input: only_pos (boolean toggle, if true, only outputs the position vector, otherwise the velocity vector is also output)
#Input: h2 (specific angular momentum, used when the orbit is parabolic, units: --> AU^2/day)
#Output: x, y, z (postion vector in the orbital plane, z will inherently always be zero, units: AU)
#Output: vx, vy, vz (velocity vector in the orbital plane, vz will inherently always be zero, units: AU/day)
def orb_to_cart_in_plane(mu, a, e, E=None, f=None, only_pos=True, h2=None):
    test_var = E if f is None else f #variable to test for length
    z = np.zeros(len(test_var)) if hasattr(test_var, '__len__') else 0 #z-axis distance is always zero in the orbital plane's reference frame (z also doubles as vz here)

    if e < 1: #elliptical orbit
        if E is None: E = compute_eccentric_anomaly(f, e)

        #intermediate values
        sinE = np.sin(E)
        cosE = np.cos(E)
        sqe = np.sqrt(1 - e*e)

        #compute orbital plane positions
        x = a*(cosE - e)
        y = a*sqe*sinE

        if only_pos: return (x, y, z)

        #intermediate values (velocity only)
        n = np.sqrt(mu / (a**3)) #mean motion
        an = a*n
        vel_denom = 1 - e*cosE

        #compute orbital plane velocities
        vx = -an*sinE / vel_denom
        vy = an*sqe*cosE / vel_denom
    elif e > 1: #hyperbolic orbit
        if f is None: f = compute_true_anomaly(E, e)

        sinf = np.sin(f)
        cosf = np.cos(f)

        p = a*(1-e*e) #semi-latus rectum
        denom = 1 + e*cosf
        
        r = p / denom
        x = r*cosf
        y = r*sinf

        if only_pos: return (x, y, z)

        vel_angle = np.arctan(e*sinf / denom) - f
        v = np.sqrt((mu / p)*(1 + e*e + 2*e*cosf))
        vx = v*np.sin(vel_angle)
        vy = v*np.cos(vel_angle)
    else: #parabolic orbit
        if f is None: f = compute_true_anomaly(E, e)
        else: #ensure f is within the range of -pi to pi (if not done, the results will be incorrect)
            if hasattr(f, '__len__'): f[f > np.pi] -= 2*np.pi #array
            else: f -= 2*np.pi * (f > np.pi) #single value

        r = h2 / (mu * (1 + cosf))
        x = r*cosf
        y = r*sinf

        if only_pos: return (x, y, z)

        v = np.sqrt(2*mu / r)
        vel_angle = f / 2 #reason for why f requires a strict range
        vx = -v*np.sin(vel_angle)
        vy = v*np.cos(vel_angle)

    return (x, y, z, vx, vy, z)
    
#Converts the position and (optionally) velocity in the orbital plane to the heliocentric ecliptic frame via a rotation matrix (from assignment 5)
#Input: i (inclination, units: radians)
#Input: b_omega (longitude of the ascending node, units: radians)
#Input: s_omega (argument of perihelion, units: radians)
#Input: x, y, z (postion vector in the orbital plane, z will inherently always be zero, units: AU)
#Input: vx, vy, vz (velocity vector in the orbital plane, vz will inherently always be zero, units: AU/day)
#Output: x_new, y_new, z_new (postion vector in the heliocentric ecliptic frame, units: AU)
#Output: vx_new, vy_new, vz_new (velocity vector in the heliocentric ecliptic frame, units: AU/day)
def orb_plane_to_helio(i, b_omega, s_omega, x, y, z, vx=None, vy=None, vz=None):
    si = np.sin(i)
    sb = np.sin(b_omega)
    ss = np.sin(s_omega)
    
    ci = np.cos(i)
    cb = np.cos(b_omega)
    cs = np.cos(s_omega)

    cbcs = cb*cs
    sbcs = sb*cs
    cbss = cb*ss
    sbss = sb*ss

    #rotation first occurs about the z-axis by s_omega, then the x-axis by i, then the z-axis again by b_omega (orbital plane axes)
    rot_m = np.array([[cbcs - sbss*ci , -cbss - sbcs*ci, sb*si], [sbcs + cbss*ci, -sbss + cbcs*ci, -cb*si], [si*ss, si*cs, ci]])
    
    pos_init = np.vstack((x, y, z)) if hasattr(x, '__len__') else np.array([x, y, z]) #handles both arrays and single values of x, y, and z
    x_new, y_new, z_new = np.matmul(rot_m, pos_init) #perform the rotation of postions

    if vx is None: return (x_new, y_new, z_new) #only return the heliocentric ecliptic postion
    
    vel_init = np.vstack((vx, vy, vz)) if hasattr(x, '__len__') else np.array([vx, vy, vz]) #handles both arrays and single values of vx, vy, and vz
    vx_new, vy_new, vz_new = np.matmul(rot_m, vel_init) #perform the rotation of velocities

    return (x_new, y_new, z_new, vx_new, vy_new, vz_new) #return the heliocentric ecliptic postion and velocity

#Function for solving Kepler's equations using a binary search approach (weird but since both functions are monotomic, it works)
#Input: e (eccentricity)
#Input: M (mean anomaly)
#Output: E (eccentric anomaly)
def _solve_kepler_eq(e, M):
    eps_lim = 10*max(float_info.epsilon, np.abs(M)*float_info.epsilon)

    if e < 1: #elliptic orbit
        if e == 0: return M #trivial case
        kepler_func = lambda e, E: E - e*np.sin(E)
        E = M #starting guess
        E_mult = np.sqrt(2) #don't scale as fast
    else: #hyperbolic orbit only, function should not be called with e = 1
        if M == 0: return 0 #trivial case
        kepler_func = lambda e, E: e*np.sinh(E) - E
        E = np.sign(M) #starting guess
        E_mult = 2
        
    min_bound = 0
    max_bound = 0

    #initialize min and max bounds
    if M < 0:
        while True:
            M_test = kepler_func(e, E)
            M_diff = M - M_test
            if np.abs(M_diff) < eps_lim: return E
            if M_diff > 0:
                min_bound = E
                break
            else: 
                max_bound = E
                E *= E_mult
    else:
        while True:
            M_test = kepler_func(e, E)
            M_diff = M - M_test
            if np.abs(M_diff) < eps_lim: return E
            if M_diff > 0:
                min_bound = E
                E *= E_mult
            else:
                max_bound = E
                break

    while True: #perform a binary search to solve for E
        E = (max_bound + min_bound) / 2 #take E at the midpoint
        M_test = kepler_func(e, E)
        M_diff = M - M_test
        if np.abs(M_diff) < eps_lim: return E
        if M_diff > 0: min_bound = E
        else: max_bound = E

#Helper function for solving Kepler's equation to allow for both single variable and array inputs for M
#See the function above for input and output variable definitions
def solve_kepler_eq(e, M):
    if hasattr(M, '__len__'):
        E = np.empty(len(M))
        for i, ma in enumerate(M): E[i] = _solve_kepler_eq(e, ma)
    else: E = _solve_kepler_eq(e, M)
    return E

#Computes the Sun-object vector at a given Julian Date
#Input: object_params (orbital elements of the target object (planet, asteroid, etc.), including the mean anomaly at the provided epoch)
#Input: epoch (the Julian Date associated with the provided orbital elements)
#Input: JD (the Julian Date at which to compute the object's position in the heliocentric ecliptic frame)
#Input: mu (the standard gravitational parameter for the object-Sun system, or just that of the Sun)
#Input: h2 (specific angular momentum, used when the orbit is parabolic, units: --> AU^2/day)
#Output: the Sun-object vector in the heliocentric ecliptic frame
def comp_sun_obj_vec(object_params, epoch, JD, mu, h2=None):
    mean_anomaly = (JD - epoch) * np.sqrt(mu / np.abs(object_params[0])**3) + object_params[-1] #mean anomaly of observation --> M = (t - t0) * n + M0

    if object_params[1] != 1: #elliptic or hyperbolic orbit
        E = solve_kepler_eq(object_params[1], mean_anomaly) #solve for E via Kepler's equation
        x, y, z = orb_plane_to_helio(*object_params[2:-1], *orb_to_cart_in_plane(mu, *object_params[:2], E=E)) #position of Earth in the heliocentric ecliptic frame
    else: #parabolic orbit
        #Barker's equation from https://en.wikipedia.org/wiki/Parabolic_trajectory (original source is the same as the one mentioned at the top of the file)
        A = 3/2 * mean_anomaly
        B = (A + np.sqrt(A*A + 1))**(1/3)
        f = 2*np.arctan(B - 1/B)
        x, y, z = orb_plane_to_helio(*object_params[2:-1], *orb_to_cart_in_plane(mu, *object_params[:2], f=f, h2=h2))
    
    return np.vstack((x, y, z)) if hasattr(x, '__len__') else np.array([x, y, z]) #return the Sun-object vector(s)

#Computes the geocentric and heliocentric distances and distance vectors in the ecliptic frame for the 3 observations, starting with the
#observational unit vectors and their corresponding Earth-Sun vectors, along with approximations of the sector ratios
#distance calculations are performed in the Cunningham frame with all geocentric and heliocentric ecliptic vectors transformed
#into and out of the Cunningham frame via rotation matrices
#Input: pts --> p1, p2, p3 (observational unit vectors in the geocentric ecliptic frame, input as a matrix)
#Input: pts_ES (3 Earth-Sun vectors in the heliocentric ecliptic frame, input as a matrix --> units: AU)
#Input: a1, a3 (sector ratio approximations --> last (a1) and first (a3) triangular area over full triangular area)
#Output: dists_geo, dists_helio (the distances to the object for each observation in the geocentric and heliocentric ecliptic frame --> units: AU)
#Output: vectors_geo, vectors_helio (the distance vectors to the object for each observation in the geocentric and heliocentric ecliptic frame --> units: AU)
#Output: nu_2 (pts_ch[2, 1] --> the z component of the unit vector for the second observation in the Cunningham frame)
def comp_dists(pts, pts_ES, a1, a3):
    err_msg = "Distances cannot be computed with this triplet of observations as "

    eta = np.cross(pts[:, 0], np.cross(pts[:, 2], pts[:, 0])) #compute the eta vector
    mag_eta = np.sqrt(eta.dot(eta))

    if mag_eta == 0: #this should never happen
        print(err_msg + "eta = 0. Are you sure that isn't a star?")
        return (None,)*5

    y_eta = eta / mag_eta #normalize the eta vector
    z_zeta = np.cross(pts[:, 0], y_eta) #p1 is also the xi unit vector (x-axis) | p1 and eta are perpendicular, therefore zeta is already normalized

    ch_frame_rot = np.vstack((pts[:, 0], y_eta, z_zeta)) #the rotation matrix for the Cunningham frame

    pts_ch = np.matmul(ch_frame_rot, pts) #perform the coord transformation for all points at once
    pts_ESC = np.matmul(ch_frame_rot, pts_ES) #convert Earth-Sun vectors to the Cunningham frame

    #catch divide by zero errors before they happen
    if pts_ch[2, 1] == 0: #this might happen rarely
        print(err_msg + "nu_2 = 0.")
        return (None,)*5
    if a1 == 0 or a3 == 0: #this really should never happen
        if a1 == 0 and a3 == 0: amsg = "a1 = 0 and a3 = 0"
        elif a1 == 0: amsg = "a1 = 0"
        else: amsg = "a3 = 0"
        print(err_msg + amsg + ". Ideally observations should be taken at different times.")
        return (None,)*5

    #compute the distances to the objects from the Earth in the Cunningham frame
    p1, p3, p2 = np.dot(np.array([a1, -1, a3]), pts_ESC.T)
    p2 /= -pts_ch[2, 1]
    p3 = (p3 + p2*pts_ch[1, 1]) / (a3*pts_ch[1, 2])
    p1 = (p1 + p2*pts_ch[0, 1] - a3*p3*pts_ch[0, 2]) / a1

    #a rotation of 180 degrees for all three observations produces the same vectors_geo as with no rotation, due to both pts and dists_geo being 
    #effectively multiplied by negative one, thefore if p1 is negative, all distances need to be flipped from negative to positive
    dists_geo = np.array([p1, p2, p3])*np.sign(p1)
    vectors_geo = dists_geo*pts
    vectors_helio = vectors_geo - pts_ES
    dists_helio = np.linalg.norm(vectors_helio, axis=0)

    return (dists_geo, dists_helio, vectors_geo, vectors_helio, pts_ch[2, 1]) #distances and vectors in the geocentric and heliocentric ecliptic frame (and nu_2)

#Compute the heliocentric (ecliptic frame) velocities of the object at each of the 3 observational points
#Input: dists_helio (distances to the object for each observation --> units: AU)
#Input: vectors_helio (heliocentric postion vectors for the object --> units: AU)
#Input: mu (standard gravitational parameter of the Sun, object is assumed to be much less massive than the Sun)
#Input: tau1, tau2, tau3 (time deltas for each pair of observations --> units: days)
#Ouput: vels_helio (heliocentric velocity vectors --> AU/day)
def comp_vels_helio(dists_helio, vectors_helio, mu, tau1, tau2, tau3):
    taus = np.array([tau3, tau1, -tau2]) #ensures all 3 pairs of points are used to find the velocities --> (t2,t1), (t3,t2), (t1,t3)
    sigma = mu / dists_helio**3
    taus_sq = taus*taus
    f = 1 - 0.5*sigma*taus_sq
    g = taus - sigma * taus*taus_sq / 6
    return (vectors_helio[:, [1,2,0]] - f*vectors_helio) / g #vels_helio (1, 2, 0 matches the first time values from a few lines above --> t2, t3, t1)

#Computes the covariance matrix for the orbital elements (for f, only the samples for the middle observation are used, otherwise all samples are used)
#Optionally computes the mean and standard deviation of each orbital element
#Input: orbital_elems (2d array of orbital elements, one row is one sample, each observation has multiple samples)
#Input: mid_obs_inxs (the indices (typically) corresponding to the middle observation)
#Input: elems_inx (the number of samples of the orbital elements)
#Input: incl_f (boolean flag to include the true anomaly in the covariance matrix)
#Input: out_mean, out_std (boolean flags for outputing the mean and standard deviation of each orbital element as well)
#Output: cov_matrix (the covariant matrix of the orbital elements)
#Output: means (the mean for each orbital element)
#Output: standard_devs (the standard deviation for each orbital element)
def comp_cov_matrix(orbital_elems, mid_obs_inxs, elems_inx, incl_f=True, out_mean=True, out_std=True):
    cov_size = 5 + incl_f
    cov_matrix = np.empty([cov_size, cov_size])

    means_all_obs = np.mean(orbital_elems[:elems_inx, :5], axis=0) #mean values of the orbital elements excluding f, using all observations
    all_elems_sub_mean = orbital_elems[:elems_inx, :5] - means_all_obs

    cov_matrix[:5, :5] = np.matmul(all_elems_sub_mean.T, all_elems_sub_mean) / (elems_inx - 1) #may seem inefficient, but numpy optimizes symmetric matrix multiplication with a special BLAS function call, neat

    if incl_f:
        means_mid_obs = np.mean(orbital_elems[mid_obs_inxs], axis=0) #mean values of all orbital elements, using only the middle observation (where f should be the same)
        mid_elems_sub_mean = orbital_elems[mid_obs_inxs] - means_mid_obs

        cov_matrix[-1] = np.dot(mid_elems_sub_mean.T, mid_elems_sub_mean[:, -1]) / (len(mid_obs_inxs) - 1) #fill the last row with covariance (and variance) data for f (true anomaly)
        cov_matrix[:-1, -1] = np.copy(cov_matrix[-1, :-1]) #copy the last row into the last column (bottom right entry of the matrix is not copied)

    if not (out_mean or out_std): return cov_matrix

    if out_mean:
        means = np.empty(6)
        means[:-1] = means_all_obs
        means[-1] = means_mid_obs[-1] if incl_f else orbital_elems[mid_obs_inxs][0,-1] #when include f is False, the mean of f is just the value of f as there is only 1 sample

    if out_std:
        standard_devs = np.zeros(6)
        standard_devs[:cov_size] = np.sqrt(np.diag(cov_matrix)) #when include f is False, the standard deviation for f is simply zero

    if not out_std and out_mean: return (cov_matrix, means)
    if not out_mean and out_std: return (cov_matrix, standard_devs)

    return (cov_matrix, means, standard_devs) #default option

#Return the sin and cosine of Earth's tilt
@functools.lru_cache(maxsize=1)
def get_earth_tilt_trig():
    earth_tilt = np.deg2rad(23.439291111111) #Earth's axial tilt (epsilon)
    return (np.cos(earth_tilt), np.sin(earth_tilt))

#Converts arrays of Right Ascension and Declination observational values into a matrix of x, y, z geocentric ecliptic unit vectors
#Input: ra (Right Ascension --> units: degrees)
#Input: dec (Declination --> units: degrees)
#Output: x, y_ec, z_ec (geocentric ecliptic unit vector components)
def comp_obs_uvec(ra, dec):
    #convert from degrees to radians
    ra = np.deg2rad(ra)
    dec = np.deg2rad(dec)

    r = np.cos(dec) #projection of the unit observation vector onto the x-y geocentric equatorial plane
    
    #geocentric equatorial unit vector
    x = r*np.cos(ra)
    y = r*np.sin(ra)
    z = np.sin(dec)

    cos_eps, sin_eps = get_earth_tilt_trig()

    #convert the equatorial plane vector to a geocentric ecliptic vector via a simplified rotation matrix (x doesn't change)
    y_ec = y*cos_eps + z*sin_eps
    z_ec = z*cos_eps - y*sin_eps

    if hasattr(x, '__len__'): return np.vstack((x, y_ec, z_ec))
    else: return np.array([x, y_ec, z_ec])

#Computes approximations to the sector ratios (the triangular ratios) based on the differences in observational times
#Input: t1, t2, t3 (the observations times --> units: Julian Days)
#Output: a1, a3 (the approximate sector ratios, time ratios are used instead due to Kepler's second law)
def comp_taus(t1, t2, t3):
    tau1, tau2, tau3 = t3-t2, t3-t1, t2-t1
    return (tau1 / tau2, tau3 / tau2)

#Helper function to return the standard gravitational parameter of the Sun and Earth-Sun system
def get_mus():
    #mass units: solar masses, postion units: AU, velocity units: AU/day
    k = 0.01720209895 #Gauss’ gravitational constant
    G = k*k #Gauss’ gravitational constant squared
    m_Sun = 1 #assume the mass of the orbiting object is much less than that of the Sun
    m_Earth = 3.0404327497692654e-06 #mass of the Earth in Solar masses
    return G*m_Sun, G*(m_Sun+m_Earth) #standard gravitational parameters

#Function to handle the catching of edge cases for when most observations don't produce valid orbital elements (extremely rare)
#Input: samples_per_obs (array with num_obs entries, whose values are the same under normal circumstances)
#Input: elems_inx (index for accessing the orbital elements 2d array)
#Input: num_obs (number of observations)
#Output: reduction_level (flag with 3 states to determine what should be computed from the orbital elements --> e.g. covariant matrix)
#Output: obs_inx (the index of the observation to use for the true anomaly --> default if all is normal will be the middle observation)
def catch_obs_edge_cases(samples_per_obs, elems_inx, num_obs):
    max_samples_inx = np.argmax(samples_per_obs)
    max_samples = samples_per_obs[max_samples_inx]
    f_text = "The true anomaly is listed for observation "
    reduction_level = 0 #no reduction by default, output the full 6x6 covariance matrix along with the mean and standard deviation
    mid_obs_inx = num_obs // 2 #mid observation index (the observation requested for use in computing f)

    if elems_inx == 1 or max_samples == 1:
        f_text_start = "The middle observation has no valid samples. "

        if elems_inx == 1:
            print(name + " has only 1 sample for which a combination of 3 observations provides orbital elements. A covariance matrix cannot be produced.")
            reduction_level = 2 #reduced to only mean values for one sample (the mean is just the values themselves --> std of 0)
        elif max_samples == 1:
            print(name + " has at most 1 sample of the orbital elements for a given observation. The covariance matrix excludes f (the true anomaly).")
            reduction_level = 1 #reduced to a 5x5 covariance matrix (f is excluded)
    else: f_text_start = "The observation with the most valid samples is not the middle observation. "

    if samples_per_obs[mid_obs_inx] != max_samples:
        f_text = f_text_start + f_text + str(max_samples_inx+1) + "."
        obs_inx = max_samples_inx
    else:
        f_text += str(mid_obs_inx+1) + "."
        obs_inx = mid_obs_inx

    print(f_text)

    return (reduction_level, obs_inx)

#Compute the absolute or apparent magnitude and optionally the diamter of the object
#Input: mag (apparent or absolute magnitude)
#Input: albedo (assumed albedo of the object based on its orbit / object type)
#Input: mean_dists_geo (average computed geocentric distances for each obseravation)
#Input: mean_dists_helio (average computed heliocentric distances for each obseravation)
#Input: mean_vecs_geo (average computed geocentric vectors for each obseravation)
#Input: mean_vecs_helio (average computed heliocentric vectors for each obseravation)
#Input: out_abs (boolean flag to set the output to be either the absolute or the apparent magnitude)
#Output: out_mag (absolute or apparent magnitude) 
#Output: diam (diameter)
def comp_mag_diam(mag, mean_dists_geo, mean_dists_helio, mean_vecs_geo, mean_vecs_helio, out_abs=True, albedo=None):
    mean_dists_mult = mean_dists_geo*mean_dists_helio
    alpha = np.arccos(np.einsum('ij,ij->j', mean_vecs_geo, mean_vecs_helio) / mean_dists_mult)
    
    G = 0.15 #phase parameter assumed value
    A1, A2, B1, B2 = 3.33, 1.87, 0.63, 1.22 #asteroid phase function coefficients
    tan_half_alpha = np.tan(alpha/2)
    
    phi1 = np.exp(-A1*tan_half_alpha**B1)
    phi2 = np.exp(-A2*tan_half_alpha**B2)

    log_mul_dist = 5*np.log10(mean_dists_mult)
    apf = 2.5*np.log10((1-G)*phi1 + G*phi2) #asteroid phase function
    
    if out_abs: out_mag = np.mean(mag - log_mul_dist + apf) #compute the mean absolute magnitude
    else: out_mag = mag + log_mul_dist - apf #compute the apparent magnitude

    if albedo is not None:
        diam = 1329 * 10**(-out_mag / 5) / np.sqrt(albedo)
        return (out_mag, diam)
    else: return out_mag

#plot all the planets and create 3 pairs of graphs
#Input: planets (the data for all the planet (and Pluto) --> orbital elements, name, plot color)
#Input: mu (standard gravitational parameter of the Sun)
#Input: num_objs (number of unknown objects to plot later, each one gets its own plot with all the planets)
#Input: axis_label_font_size (self explanatory, needed later in the main function as well)
#Input: orbit_plot_pts (number of points (true anomalies) at which to compute the orbital position)
#Input: today (today's Julian Date)
#Input: epoch_planets (the epoch for the planet data)
#Output: all_axes (the plots, which will be added to later)
def plot_planets(planets, mu, num_objs, axis_label_font_size, orbit_plot_pts, today, epoch_planets):
    fs = np.linspace(0, 2*np.pi, orbit_plot_pts) #true anomalies at which to compute the x, y heliocentric coords
    
    num_figs = (num_objs+1) // 2
    all_axes = [None]*(2*num_figs)
    
    for i in range(num_figs):
        fig, axes = plt.subplots(1, 2) #add subplot_kw=dict(projection='3d') for 3d
        fig.subplots_adjust(left=0.05, bottom=0, right=0.95, top=1, wspace=0.15, hspace=0)
        
        all_axes[2*i] = axes[0]
        all_axes[2*i+1] = axes[1]
    
    for name, planet_data in planets.items():
        c = planet_data[0]
        planet_params = planet_data[1]

        x, y, _ = orb_plane_to_helio(*planet_params[2:-1], *orb_to_cart_in_plane(mu, *planet_params[:2], f=fs))
        xte, yte, _ = comp_sun_obj_vec(planet_params, epoch_planets, today, mu) #the positions of the planets today
        for ax in all_axes:
            ax.plot(x, y, c=c, label=name)
            ax.scatter(xte, yte, c=c, s=30, zorder=10) #plot the location of the Earth

    for ax in all_axes:
        #add a full-plot crosshair going through (0, 0) to make the circular deviation of minimally eccentric orbits easier to see
        ax.axhline(y=0, linewidth=1, color='k')
        ax.axvline(x=0, linewidth=1, color='k')

        ax.scatter(0, 0, c="yellow", s=5, zorder=10) #add a dot for the Sun (on top of the crosshairs)

        ax.set_aspect('equal', adjustable='box') #make sure the plots remain square even if the window size is adjusted
        ax.set_xlabel("X-Axis Distance from Sun (AU)", fontsize=axis_label_font_size)
        ax.set_ylabel("Y-Axis Distance from Sun (AU)", fontsize=axis_label_font_size)

    return all_axes

#Initalized the plots for Dec vs RA and Mag vs Time
#Input : axis_label_font_size (self explanatory, needed later in the main function as well)
#Input: num_objs (number of unknown objects to plot later, each one gets its own plot with all the planets)
#Output: combined_fig (the figure containing 2 plots, needed later for date formating)
#Output: separate_fig (the figure containing pairs of plots for each object, needed later for date formating)
#Output: combined_obj_axes (the plots where all object are plotted together, which will be added to later)
#Output: separate_obj_axes (the plots where all object are plotted separately, which will be added to later)
#Outpit: colors (the colors for use when plotting the objects)
def initialize_radecmag_plots(axis_label_font_size, num_objs):
    combined_fig, combined_obj_axes = plt.subplots(1, 2) #axes containing the combined plots for dec vs ra and mag vs time
    separate_fig, separate_obj_axes = plt.subplots(2, num_objs) #axes containing the separate plots for dec vs ra and mag vs time

    combined_fig.subplots_adjust(left=0.05, bottom=0.15, right=0.95, top=0.95, wspace=0.3, hspace=None)
    separate_fig.subplots_adjust(left=0.05, bottom=0.12, right=0.95, top=0.95, wspace=0.4, hspace=0.2)

    #set titles for the combined plot (separate plots won't have titles as it would be too crowded)
    combined_obj_axes[0].set_title("Paths of Observed Objects as Seen From Earth", fontsize=axis_label_font_size + 2)
    combined_obj_axes[1].set_title("Brightness over Time", fontsize=axis_label_font_size + 2)
    
    for ax in chain(separate_obj_axes[0, :], [combined_obj_axes[0]]): #Dec vs Ra
        ax.set_xlabel("Right Ascension (deg)", fontsize=axis_label_font_size)
        ax.set_ylabel("Declination (deg)", fontsize=axis_label_font_size)
    

    for ax in chain(separate_obj_axes[1, :], [combined_obj_axes[1]]): #Mag vs Julian Date
        ax.set_xlabel("Time (Julian Date)", fontsize=axis_label_font_size)
        ax.set_ylabel("Apparent Magnitude", fontsize=axis_label_font_size)

        ax.xaxis.set_major_formatter(mdates.DateFormatter('%m/%d/%Y'))
        ax.xaxis.set_major_locator(mdates.AutoDateLocator())

    #supports up to ten objects (there are better ways of doing this, but this works for now)
    colors = ["r", "g", "b", "c", "m", "k", "brown", "teal", "grey", "purple"][:num_objs]

    return (combined_fig, separate_fig, combined_obj_axes, separate_obj_axes, colors)

  
if __name__ == "__main__":
    asteroids = initialize_obs_data()
    planets, epoch_planets = initialize_planet_data()
    GMS, GMSE = get_mus()
    today = to_datetime('now').to_julian_date() #today's Julian Date

    #initialize plots
    num_objs = len(asteroids)
    axis_label_font_size = 12
    orbit_plot_pts = 200 #number of samples of f used for plotting the orbits
    orbit_axes = plot_planets(planets, GMS, num_objs, axis_label_font_size, orbit_plot_pts, today, epoch_planets) #plot the planets now, and return the axes
    combined_fig, separate_fig, combined_obj_axes, separate_obj_axes, colors = initialize_radecmag_plots(axis_label_font_size, num_objs)
    
    num_obs = len(list(asteroids.values())[0]) #number of observations (assumed the same for each object)
    pts_combs = np.array(list(combinations(np.arange(num_obs), 3))) #all possible combinations of 3 pts (if num_obs is much larger, this can be changed to a subsample)
    num_samples = 3*len(pts_combs) #(maximum) number of samples taken of the orbital elements

    #orbital element output string variables
    osf = " = {:0.3f} ({:0.3f}) " #output string format
    elem_names = ["a", "e", "i", "b_omega", "s_omega", "f"]
    elem_unit_names = ["AU", "", "deg", "deg", "deg", "deg"]

    for ii, (name, obs) in enumerate(asteroids.items()): #loop over all asteroids / comets / planets?
        pts = comp_obs_uvec(obs[:, 1], obs[:, 2])
        pts_ES = -1*comp_sun_obj_vec(planets.get("Earth")[1], epoch_planets, obs[:, 0], GMSE) #Earth-Sun vectors

        orbital_elems = np.empty([num_samples, 6]) #a, e, i, b_omega, s_omega, f (degrees)
        elems_inx = 0 #index for accessing the orbital elements 2d array (the final value equals the total number of valid orbital element sets)
        obs_inxs = [[] for _ in range(num_obs)] #collection of indices of orbital_elems corresponding to a given observation (where f (true anomaly) should be the same)
        news = 0 #average value of all nus, weighted based on number of valid orbital element sets per nu

        #mean values for the distances (and distance vectors) for use in computing the object's absolute magnitude
        mean_dists_geo = np.zeros(num_obs)
        mean_dists_helio = np.zeros(num_obs)
        mean_vecs_geo = np.zeros([3, num_obs])
        mean_vecs_helio = np.zeros([3, num_obs])

        angular_momentum = np.zeros(3) #used only by parabolic orbits (so it will essentially never be used, but just in case)

        print("\n\nObject: " + name)
        for triplet_inx, pts_inxs in enumerate(pts_combs): #iterate over all possible combinations of 3 observations
            #observational timing related info
            t1, t2, t3 = obs[pts_inxs, 0] #observation times
            tau1, tau2, tau3 = t3-t2, t3-t1, t2-t1 #observation time deltas
            
            if tau2 == 0:
                print("Distances cannot be computed with this triplet of observations as tau2 = 0.")
                continue

            a1, a3 = tau1 / tau2, tau3 / tau2 #sector ratio approximations
            
            dists_geo, dists_helio, vectors_geo, vectors_helio, nu_2 = comp_dists(pts[:, pts_inxs], pts_ES[:, pts_inxs], a1, a3)
            if dists_geo is None: continue
            vels_helio = comp_vels_helio(dists_helio, vectors_helio, GMS, tau1, tau2, tau3)
            
            print("Observation Triplet #" + str(triplet_inx+1) + " | Observations: ", [x+1 for x in pts_inxs])
            print("Geocentric Distances: P1 = {:0.3f} AU, P2 = {:0.3f} AU, P3 = {:0.3f} AU".format(*dists_geo))
            print("Heliocentric Distances: P1 = {:0.3f} AU, P2 = {:0.3f} AU, P3 = {:0.3f} AU".format(*dists_helio))
            print("nu_2 = {:0.8f}".format(nu_2))
            print("a1 = {:0.3f}, a3 = {:0.3f}".format(a1, a3))

            for inx, pts_inx in enumerate(pts_inxs):
                #requires too many logic changes to the subroutine to allow for a simple matrix input
                a, e, i, b_omega, s_omega, f = compute_orbital_elements(GMS, *vectors_helio[:, inx], *vels_helio[:, inx])
                if a is None: continue

                #only include the angular momentum in the eventual average value, if the orbital elements where valid
                angular_momentum += np.cross(vectors_helio[inx], vels_helio[inx]) #only need for parabolic orbits, but the type isn't know until the mean eccentricity is computed
                
                #only include mean distances if the orbital elements for that point could be computed
                mean_dists_geo[pts_inx] += dists_geo[inx]
                mean_dists_helio[pts_inx] += dists_helio[inx]
                mean_vecs_geo[:, pts_inx] += vectors_geo[:, inx]
                mean_vecs_helio[:, pts_inx] += vectors_helio[:, inx]

                news += nu_2

                orbital_elems[elems_inx] = np.array([a, e, i, b_omega, s_omega, f])
                obs_inxs[pts_inx].append(elems_inx) #orbital elements where f should be the same
                elems_inx += 1

            print()
        
        if elems_inx == 0: #reduction level 3 (reduced to nothing)
            print(name + " has no combinations of 3 observations which can provide orbital elements.\n")
            continue

        #compute the mean values of the distances and velocities
        samples_per_obs = np.array([len(x) for x in obs_inxs]) #should be an array with num_obs entries all equal to the same value (under normal circumstances)
        non_zero_samples = np.nonzero(samples_per_obs) #ignore observations that had no valid samples (extremely unlikely this will ever occur, but it could)
        mean_dists_geo = mean_dists_geo[non_zero_samples] / samples_per_obs[non_zero_samples]
        mean_dists_helio = mean_dists_helio[non_zero_samples] / samples_per_obs[non_zero_samples]
        mean_vecs_geo = mean_vecs_geo[:, non_zero_samples] / samples_per_obs[non_zero_samples]
        mean_vecs_helio = mean_vecs_helio[:, non_zero_samples] / samples_per_obs[non_zero_samples]

        num_valid_obs = len(non_zero_samples[0])
        mean_vecs_shape = (3, num_valid_obs)
        mean_vecs_geo.shape = mean_vecs_shape
        mean_vecs_helio.shape = mean_vecs_shape
        
        news /= elems_inx #weighted average value of nu_2

        means_dist_str = ", ".join(["O" + str(x+1) + " = {:0.3f} AU" for x in range(num_valid_obs)])
        print("Mean Geocentric Distances: " + means_dist_str.format(*mean_dists_geo))
        print("Mean Heliocentric Distances: " + means_dist_str.format(*mean_dists_helio))
        print("Mean value of nu_2: {:0.8f}".format(news))
        print()
        
        reduction_level, obs_inx = catch_obs_edge_cases(samples_per_obs, elems_inx, num_obs)

        #compute the convariance matrix, means, and standard deviations of the orbital elements
        if reduction_level == 2:
            cov_matrix = None #no covariance matrix
            means = orbital_elems[0] #only one set of orbital elements
            standard_devs = np.zeros(6) #set the standard deviations to zero
        else: cov_matrix, means, standard_devs = comp_cov_matrix(orbital_elems, obs_inxs[obs_inx], elems_inx, not reduction_level)

        if cov_matrix is not None:
            print("\nCovariance matrix: ")
            print(cov_matrix)
            print()

        mean_std = np.empty(12) #interleaved mean and standard deviation values, used for output to terminal
        mean_std[0::2] = means
        mean_std[1::2] = standard_devs

        means[2:] *= np.pi/180 #now that means isn't being output to the terminal, convert all degrees to radians

        orb_out_str = ""
        for elem_name, elem_unit_name in zip(elem_names, elem_unit_names): orb_out_str += elem_name + osf + elem_unit_name + ", "
        orb_out_str = orb_out_str[:-2]
        print(orb_out_str.format(*mean_std))
        print()

        #magnitude calculation section
        albedos = np.arange(0.02, 0.32, 0.02) #range of albedos with which possible diameters are calculated

        H_mag, diams = comp_mag_diam(obs[non_zero_samples, -1], mean_dists_geo, mean_dists_helio, mean_vecs_geo, mean_vecs_helio, albedo=albedos)
        print("Absolute magnitude: {:0.3f}".format(H_mag))
        print()
        for albedo, diam in zip(albedos, diams): print("Albedo = {:0.3f} | Diameter = {:0.3f} km".format(albedo, diam))
        print()

        #ephemeris calculation section
        obs_valid_inx = non_zero_samples[0][0] #index of the first observation with valid orbital elements (should be obs 1)
        f_first = np.deg2rad(np.mean(orbital_elems[obs_inxs[obs_valid_inx]][:, -1])) #average value of the true anomaly for the first observation with valid orbital elements (again should be obs 1)
        E = compute_eccentric_anomaly(f_first, means[1]) #convert f to E
        epoch = obs[obs_valid_inx, 0] #starting date
        JDs = np.arange(91) + epoch + 30 #array of Julian Dates for 90 days starting 30 days after the first valid observation
        obj_params = np.copy(means)
        
        #handle the separate conditions for the 3 types of orbit
        #replace the mean true anomaly at the (almost always) middle observation with the mean anomaly at the (almost always) first observation
        h2=None #specific angular momentum, only used by parabolic orbits
        if means[1] < 1: obj_params[-1] = E - means[1]*np.sin(E) #elliptical orbit
        elif means[1] > 1: obj_params[-1] = means[1]*np.sinh(E) - E #hyperbolic orbit
        else: #parabolic orbit
            angular_momentum /= elems_inx #average the angular momentum
            h2 = np.dot(angular_momentum, angular_momentum) #squared average angular momentum
            obj_params[0] = 2**(1/3)*h2/(2*GMS) #for use in the mean motion equation a is now 2^(1/3) times the periapsis
            obj_params[-1] = E + E**3/3

        #compute postition vectors
        pos_vecs_obj_helio = comp_sun_obj_vec(obj_params, epoch, JDs, GMS, h2=h2)
        pos_vecs_earth_helio = comp_sun_obj_vec(planets.get("Earth")[1], epoch_planets, JDs, GMSE)
        pos_vecs_obj_geo = pos_vecs_obj_helio - pos_vecs_earth_helio

        #rotate the geocentric ecliptic frame object vectors into the geocentric equatorial plane, then compute ra and dec
        cos_eps, sin_eps = get_earth_tilt_trig()
        ecl_to_equ_rot = np.array([[1, 0, 0], [0, cos_eps, -sin_eps], [0, sin_eps, cos_eps]])
        pos_vecs_obj_equ = np.matmul(ecl_to_equ_rot, pos_vecs_obj_geo)
        ra, dec = ra_and_dec(*pos_vecs_obj_equ)

        #compute apparent magnitude for each day
        dists_obj_geo = np.linalg.norm(pos_vecs_obj_geo, axis=0)
        dists_obj_helio = np.linalg.norm(pos_vecs_obj_helio, axis=0)
        app_mags = comp_mag_diam(H_mag, dists_obj_geo, dists_obj_helio, pos_vecs_obj_geo, pos_vecs_obj_helio, out_abs=False)

        #plotting section
        
        #plotting the object's orbits onto the existing plots
        fs, Es = None, None #one will be None, the other will have values
        max_perc = 0.98 #percentage of the way to the maximim allowed angle for hyperbolic and parabolic orbits
        extra_scale = 1.2 #additional amount of the Solar System to show vs some baseline AU distance
        large_scale_reduce = 10.5/12 #reduce the amount of the Solar System to show for large orbits
        min_large_AU = 470 #minimum distance to show on the plot for hyperbolic and parabolic orbits
        if means[1] < 1:
            Es = np.linspace(0, 2*np.pi, orbit_plot_pts) #elliptic orbits have their entire orbit plotted
            au_lim = means[0]*(1+means[1])*extra_scale
        elif means[1] > 1:
            f_lim = np.arccos(-1/means[1]) #maximum defined angle of f for a hyperbolic orbit
            fs = np.linspace(-f_lim*max_perc, f_lim*max_perc, orbit_plot_pts) #plots for hyperbolic orbits show most of the region where the orbit is defined
            au_lim = max(min_large_AU, means[0]*(1-means[1])*extra_scale) #scale of the graph is min_large_AU AU or extra_scale times the periapsis, whichever is largest
        else:
            fs = np.linspace(-np.pi*max_perc, np.pi*max_perc, orbit_plot_pts) #plots for parabolic orbits have their orbit cut off when f is close to apoapsis
            au_lim = max(min_large_AU, h2/(2*GMS)*extra_scale) #scale of the graph is min_large_AU AU or extra_scale times the periapsis, whichever is largest
        if au_lim > min_large_AU: au_lim *= large_scale_reduce

        custom_color = ("g" if name == "X0613" else "r")
        custom_color_inv = ("r" if name == "X0613" else "g")
        x, y, _ = orb_plane_to_helio(*means[2:-1], *orb_to_cart_in_plane(GMS, *means[:2], E=Es, f=fs))
        orbit_axes[ii].plot(x, y, c=custom_color, label=name+" Orbit") #the object assigned to me is colored green, the others are red
        orbit_axes[ii].set_title("Target: " + name, fontsize=axis_label_font_size + 2)
        orbit_axes[ii].set_xlim([-au_lim, au_lim])
        orbit_axes[ii].set_ylim([-au_lim, au_lim])

        #plot the object on the date of the first (valid) observation and on today's date
        for jd, clr in zip([epoch, today], [custom_color, custom_color_inv]):
            date_str = to_datetime(jd, origin="julian", unit='D').strftime('%m/%d/%Y')
            x, y, _ = comp_sun_obj_vec(obj_params, epoch, jd, GMS, h2=h2)
            orbit_axes[ii].scatter(x, y, c=clr, s=30, zorder=10, label=name + " on " + date_str) #plot the location of the object 
        
        #plotting dec vs ra and mag vs time
        JD_plotting = [to_datetime(JD, origin="julian", unit='D').date() for JD in JDs]
        linewidth = 4 if name == "X0613" else 2 #thicker line given to my assigned object
        
        #handle the case where RA crosses from 0 to 360 or 360 to 0 degrees (prevents unwanted lines from being drawn from one side of the plot to the other)
        discont = np.where(np.abs(np.diff(ra)) > 300)[0] #find the crossings
        start_j = 0
        for j in discont: #loop over all crossings
            combined_obj_axes[0].plot(ra[start_j:j+1], dec[start_j:j+1], c=colors[ii], linewidth=linewidth)
            separate_obj_axes[0, ii].plot(ra[start_j:j+1], dec[start_j:j+1], c=colors[ii], linewidth=linewidth)
            start_j = j+1

        combined_obj_axes[0].plot(ra[start_j:], dec[start_j:], c=colors[ii], label=name, linewidth=linewidth)
        combined_obj_axes[1].plot(JD_plotting, app_mags, c=colors[ii], label=name, linewidth=linewidth)
        separate_obj_axes[0, ii].plot(ra[start_j:], dec[start_j:], c=colors[ii], linewidth=linewidth)
        separate_obj_axes[1, ii].plot(JD_plotting, app_mags, c=colors[ii], linewidth=linewidth)
        separate_obj_axes[0, ii].set_title(name, fontsize=axis_label_font_size + 2)
    
    for ax in chain(orbit_axes, combined_obj_axes): ax.legend(loc="upper right") #add legends
    for ax in chain(separate_obj_axes[1, :], [combined_obj_axes[1]]): plt.setp(ax.xaxis.get_majorticklabels(), rotation=70) #display the dates for the mag vs time plots nicely

    plt.show()
