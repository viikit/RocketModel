import matplotlib.pyplot as plt
import numpy as np
# from scipy.integrate import solve_ivp
import scipy.integrate as sci
import time
import plotly.graph_objects as go
from scipy.interpolate import interp1d
import inputNovaRocket as input

# Input values
# constants
G = input.G

# Planet Earth
Rplanet = input.Rplanet
Mplanet = input.Mplanet
# Wind
alpha = input.alpha

# Parameters of rocket
#initial conditions for single stage rocket
x0 = input.x0
z0 = input.z0
velx0 = input.velx0
velz0 = input.velz0
period = input.period
period1 = input.period1

# Mass
m1all = input.m1all
m1prop = input.m1prop
m2all = input.m2all
m2prop = input.m2prop
m3all = input.m3all
m0 = input.m0
mfairing = input.mfairing

#thrust
thrustF1 = input.thrust1
thrust2 = input.thrust2

# Specific Impulse
Isp1 = input.Isp1
Isp2 = input.Isp2

# Duration Times in seconds
t0 = input.t0
t_vert = input.t_vert
t_pitchover = input.t_pitchover
delta_t1a = input.delta_t1a
delta_t1b = input.delta_t1b
delta_t1off = input.delta_t1off
t1end = input.t1end
delta_t02 = input.delta_t02
t2start = input.t2start
t2man = input.t2man
delta_t2 = input.delta_t2
delta_t20ff = input.delta_t20ff
t2end = input.t2end

psi1 = input.psi1
dpsi1 = input.dpsi1

psi2 = input.psi2
dpsi2 = input.dpsi2

psi3 = input.psi3
dpsi3 = input.dpsi3

psi4 = input.psi4
dpsi4 = input.dpsi4

psi5 = input.psi5
dpsi5 = input.dpsi5


mdotF1 = input.mdotF1
mdot2 = input.mdot2

m_prelaunch = input.m_prelaunch
m_launch = input.m_launch
mdry1 = input.mdry1
thrustfactor1 = input.thrustfactor1
thrustfactor2 = input.thrustfactor2
thrustfactor3 = input.thrustfactor3
thrustfactor4 = input.thrustfactor4
thrustfactor5 = input.thrustfactor5
thrustfactor6 = input.thrustfactor6

D = input.D
S = input.S
L1 = input.L1
L2 = input.L2
CD = input.CD
cdwind = input.cdwind

# WIND and Density
windfactor = input.windfactor
windspeed_file = np.loadtxt('wind-alt.txt')
alt_data = windspeed_file[:, 0]
wind_data = windspeed_file[:, 1]
wind_interp = interp1d(alt_data, wind_data, 'cubic')

density_file = np.loadtxt('density-alt.txt')
alt_dendata = density_file[:,0]
rho_dendata = density_file[:,1]
density_interp = interp1d(alt_dendata, rho_dendata, 'cubic')

# gravitational acceleration
def gravity(x, z):
    global Rplanet, Mplanet, G

    r = np.sqrt(x**2 + z**2)

    if r < Rplanet:
        accelx = 0.0
        accelz = 0.0
    else:
        accelx = G * Mplanet / (r ** 3) * x
        accelz = G * Mplanet / (r ** 3) * z

    return np.asarray([accelx, accelz]),r

def propulsion(x,z,velx,velz,t):

    r_vector = np.sqrt(x ** 2 + z ** 2)
    rvector_mag = abs(r_vector)

    Vel_mag = abs(np.sqrt(velx ** 2 + velz ** 2))

################ 0

    if t <= t_vert:
        theta = 0.0
        thrustF = thrustF1
        mdot = mdotF1
        thrustx = thrustF
        thrustz = 0.0

################ 1 pitchover

    if t > t_vert and t <= t_pitchover:
        theta = (psi1*np.pi/180) + (dpsi1*np.pi/180)*((t-t_vert) / (t_pitchover - t_vert))

        thrustF = thrustF1 * thrustfactor1

        mdot = mdotF1

        # RADIUS-VECTOR
        thrustx = thrustF * ((np.cos(theta) * x / rvector_mag) - (np.sin(theta) * z / rvector_mag))
        thrustz = thrustF * ((np.sin(theta) * x / rvector_mag) + (np.cos(theta) * z / rvector_mag))

        # VELOCITY VECTOR
        # thrustx = thrustF * ((np.cos(theta) * velx / Vel_mag) - (np.sin(theta) * velz / Vel_mag))
        # thrustz = thrustF * ((np.sin(theta) * velx / Vel_mag) + (np.cos(theta) * velz / Vel_mag))

    ################ 2 stage 1a

    if t > t_pitchover and t <= delta_t1a:
        theta = (psi2*np.pi/180) + (dpsi2*np.pi/180)*((t - t_pitchover) / (delta_t1a - t_pitchover))

        thrustF = thrustF1 * thrustfactor2

        mdot = mdotF1

        # RADIUS-VECTOR
        # thrustx = thrustF * ((np.cos(theta) * x / rvector_mag) - (np.sin(theta) * z / rvector_mag))
        # thrustz = thrustF * ((np.sin(theta) * x / rvector_mag) + (np.cos(theta) * z / rvector_mag))

        # VELOCITY VECTOR
        thrustx = thrustF * ((np.cos(theta) * velx / Vel_mag) - (np.sin(theta) * velz / Vel_mag))
        thrustz = thrustF * ((np.sin(theta) * velx / Vel_mag) + (np.cos(theta) * velz / Vel_mag))

    ################ 3 stage 1b

    if t > delta_t1a and t <= (delta_t1a + delta_t1b):
        theta = (psi3*np.pi/180) + (dpsi3*np.pi/180)*((t - delta_t1a) / delta_t1b)

        thrustF = thrustF1 * thrustfactor3

        mdot = mdotF1

        # RADIUS-VECTOR
        # thrustx = thrustF * ((np.cos(theta) * x / rvector_mag) - (np.sin(theta) * z / rvector_mag))
        # thrustz = thrustF * ((np.sin(theta) * x / rvector_mag) + (np.cos(theta) * z / rvector_mag))

        # VELOCITY VECTOR
        thrustx = thrustF * ((np.cos(theta) * velx / Vel_mag) - (np.sin(theta) * velz / Vel_mag))
        thrustz = thrustF * ((np.sin(theta) * velx / Vel_mag) + (np.cos(theta) * velz / Vel_mag))

    if t > (delta_t1a + delta_t1b) and t <= t1end:
        theta = 0.0
        thrustF = 0.0
        mdot = - (m1all + m_prelaunch + delta_t1a * mdotF1 + delta_t1b * mdotF1 + mfairing) / delta_t1off

        # RADIUS-VECTOR
        # thrustx = thrustF * ((np.cos(theta) * x / rvector_mag) - (np.sin(theta) * z / rvector_mag))
        # thrustz = thrustF * ((np.sin(theta) * x / rvector_mag) + (np.cos(theta) * z / rvector_mag))

        # VELOCITY VECTOR
        thrustx = thrustF * ((np.cos(theta) * velx / Vel_mag) - (np.sin(theta) * velz / Vel_mag))
        thrustz = thrustF * ((np.sin(theta) * velx / Vel_mag) + (np.cos(theta) * velz / Vel_mag))

    if t > t1end and t <= t2start:
        theta = 0.0
        thrustF = 0.0
        mdot = 0.0

        # RADIUS-VECTOR
        # thrustx = thrustF * ((np.cos(theta) * x / rvector_mag) - (np.sin(theta) * z / rvector_mag))
        # thrustz = thrustF * ((np.sin(theta) * x / rvector_mag) + (np.cos(theta) * z / rvector_mag))

        # VELOCITY VECTOR
        thrustx = thrustF * ((np.cos(theta) * velx / Vel_mag) - (np.sin(theta) * velz / Vel_mag))
        thrustz = thrustF * ((np.sin(theta) * velx / Vel_mag) + (np.cos(theta) * velz / Vel_mag))

    ################ 4 stage 2a - maneuver

    if t > t2start and t <= t2man:
        theta = (psi4*np.pi/180)+(dpsi4*np.pi/180)*((t-t2start)/(t2man-t2start))
        thrustF = thrust2 * thrustfactor4
        mdot = mdot2

        # RADIUS-VECTOR
        thrustx = thrustF * ((np.cos(theta) * x / rvector_mag) - (np.sin(theta) * z / rvector_mag))
        thrustz = thrustF * ((np.sin(theta) * x / rvector_mag) + (np.cos(theta) * z / rvector_mag))

        # VELOCITY VECTOR
        # thrustx = thrustF * ((np.cos(theta) * velx / Vel_mag) - (np.sin(theta) * velz / Vel_mag))
        # thrustz = thrustF * ((np.sin(theta) * velx / Vel_mag) + (np.cos(theta) * velz / Vel_mag))

    ################ 5 stage 2

    if t > t2man and t <= t2end:
        theta = (psi5*np.pi/180)+(dpsi5*np.pi/180)*((t-t2man)/(t2end-t2man))
        thrustF = thrust2 * thrustfactor5
        mdot = mdot2

        # # RADIUS-VECTOR
        # thrustx = thrustF * ((np.cos(theta) * x / rvector_mag) - (np.sin(theta) * z / rvector_mag))
        # thrustz = thrustF * ((np.sin(theta) * x / rvector_mag) + (np.cos(theta) * z / rvector_mag))

        # VELOCITY VECTOR
        thrustx = thrustF * ((np.cos(theta) * velx / Vel_mag) - (np.sin(theta) * velz / Vel_mag))
        thrustz = thrustF * ((np.sin(theta) * velx / Vel_mag) + (np.cos(theta) * velz / Vel_mag))

    if t > t2end:
        theta = 0.0
        thrustF = 0.0
        mdot = 0.0

        # RADIUS-VECTOR
        # thrustx = thrustF * ((np.cos(theta) * x / rvector_mag) - (np.sin(theta) * z / rvector_mag))
        # thrustz = thrustF * ((np.sin(theta) * x / rvector_mag) + (np.cos(theta) * z / rvector_mag))

        # VELOCITY VECTOR
        thrustx = thrustF * ((np.cos(theta) * velx / Vel_mag) - (np.sin(theta) * velz / Vel_mag))
        thrustz = thrustF * ((np.sin(theta) * velx / Vel_mag) + (np.cos(theta) * velz / Vel_mag))

## calculate thrust components


    return np.asarray([thrustx, thrustz]), mdot


def windf(x, z):
    r = np.sqrt(x ** 2 + z ** 2)
    alt = r - Rplanet
    wind = wind_interp(alt)
    wx = ((np.cos(alpha) * x / abs(r)) - (np.sin(alpha) * z / abs(r))) * wind * windfactor
    wz = ((np.sin(alpha) * x / abs(r)) + (np.cos(alpha) * z / abs(r))) * wind * windfactor
    return np.asarray([wx,wz])


# Equations of motion
# F = m*a = m * zddot
# z iz the altitude from the center of the planet alont the north pole
# x is the altitude from center through the equator
# zdot is the velocity along z
# zddot is the acceleration along z
# Second order differential equation

def Derivatives(state, t):
    # state vector
    x = state[0]
    z = state[1]
    velx = state[2]
    velz = state[3]
    m_rocket = state[4]
    # compute zdot- kinematic relationship
    zdot = velz
    xdot = velx

    # compute the total forces
    # Gravity
    accel, r = gravity(x,z)
    gravityF = - accel * m_rocket
    # Aerodynamics
    altitude = r - Rplanet
    rho0 = 1.293  # air density at sea level [kg/m3]
    R = 8.314  # universal gas constant of air [J/mol.K]
    M = 0.02897  # kg/mol
    T0 = 288.15  # K
    L = 0.0065  # K/m
    g0 = 9.81
    H0 = R * T0 / (g0 * M) - L / T0  # meters
    rho = rho0 * np.exp(-(max(1, altitude)) / H0)

    V = np.sqrt(velx ** 2 + velz ** 2)
    aeroF = - 0.5 * min(rho, 1.293) * S * abs(V) * CD * np.asarray([velx, velz])


    # Thrust
    thrustF, mdot = propulsion(x, z, velx, velz, t)

    if r < Rplanet:
        Forces = [0.0, 0.0, 0.0]
    else:
        Forces = gravityF + thrustF + aeroF

    gravity_Mag = np.sqrt(gravityF[0]**2 + gravityF[1]**2)
    thrust_Mag = np.sqrt(thrustF[0]**2 + thrustF[1]**2)
    drag_Mag = np.sqrt(aeroF[0]**2 + aeroF[1]**2)

    # compute the acceleration
    if m_rocket > 0:
        Vdot = Forces / m_rocket
    else:
        Vdot = [0.0,0.0]
        mdot = 0.0

    # compute the statedot
    statedot = np.asarray([xdot, zdot, Vdot[0], Vdot[1], mdot])
    return statedot


def Derivatives1(state, t):
    # state vector
    xstage1 = state[0]
    zstage1 = state[1]
    velxstage1 = state[2]
    velzstage1 = state[3]
    m_stage1 = 950.0
    # compute zdot- kinematic relationship
    zdot = velzstage1
    xdot = velxstage1


    # compute the total forces
    # Gravity
    accel, r = gravity(xstage1,zstage1)
    gravityF = - accel * m_stage1

    # Aerodynamics
    altitude = r - Rplanet
    rho = density_interp(altitude)
    vwind = windf(xstage1, zstage1)
    V = np.sqrt((velxstage1 - vwind[0]) ** 2 + (velzstage1 - vwind[1]) ** 2)

    aeroF = - 0.5 * rho * L1 * abs(V) * cdwind * np.asarray([velxstage1-vwind[0], velzstage1-vwind[1]])


    if r <= Rplanet:
        statedot1 = np.asarray([0.0, 0.0, 0.0, 0.0])
        return statedot1
    else:
        Forces = gravityF + aeroF #+ windForce
        # compute the acceleration
        Vdot = np.asarray(Forces) / m_stage1

        # compute the statedot1
        statedot1 = np.asarray([xdot, zdot, Vdot[0], Vdot[1]])
        return statedot1

def Derivatives1nowind(state, t):
    # state vector
    x = state[0]
    z = state[1]
    velx = state[2]
    velz = state[3]
    m_stage1 = 950.0
    # compute zdot- kinematic relationship
    zdot = velz
    xdot = velx


    # compute the total forces
    # Gravity
    accel, r = gravity(x,z)
    gravityF = - accel * m_stage1

    # Aerodynamics
    altitude = r - Rplanet
    rho = density_interp(altitude)

    V = np.sqrt(velx ** 2 + velz ** 2)
    aeroF = - 0.5 * min(rho, 1.293) * L1 * abs(V) * cdwind * np.asarray([velx, velz])

    if r <= Rplanet:
        Forces = [0.0, 0.0]
        xdot = 0.0
        zdot = 0.0
    else:
        Forces = gravityF + aeroF

    # compute the acceleration
    Vdot = np.asarray(Forces) / m_stage1

    # compute the statedot1
    statedot1nowind = np.asarray([xdot, zdot, Vdot[0], Vdot[1]])
    return statedot1nowind


# Populate initial condition vector
stateinitial = np.asarray([x0, z0, velx0, velz0, m_launch])


# Time window for the whole rocket, first stage and second stage
tout = np.linspace(0, period, 100000)
tout1 = np.linspace(0, period1, 100000)
tout1fall = np.linspace(0, t1end, 10000)


# Stateout for the rocket, first stage and second stage
stateout = sci.odeint(Derivatives, stateinitial, t=tout)
stateout1init = sci.odeint(Derivatives, stateinitial, t=tout1fall)


# first stage setup initial parameters
x1init = stateout1init[:,0]
z1init = stateout1init[:,1]
velx1init = stateout1init[:,2]
velz1init = stateout1init[:,3]


# first stage setup initial parameters
x10 = x1init[-1]
z10 = z1init[-1]
velx10 = velx1init[-1]
velz10 = velz1init[-1]

stateinitial1 = np.asarray([x10, z10, velx10, velz10])


# Stateout rocket at 1st stage separation


stateout1 = sci.odeint(Derivatives1, stateinitial1, t=tout1)
stateout1nowind = sci.odeint(Derivatives1nowind, stateinitial1, t=tout1)


# Rename variables
xout = stateout[:, 0]
zout = stateout[:, 1]
altitude = np.sqrt(xout**2 + zout**2) - Rplanet
velxout = stateout[:, 2]
velzout = stateout[:, 3]
velout = np.sqrt(velxout**2 + velzout**2)
massout = stateout[:, 4]

# Rename variables 1st stage
xout1 = stateout1[:, 0]
zout1 = stateout1[:, 1]
altitude1 = np.sqrt(xout1**2 + zout1**2) - Rplanet
velxout1 = stateout1[:, 2]
velzout1 = stateout1[:, 3]
velout1 = np.sqrt(velxout1**2 + velzout1**2)


# Rename variables 1st stage no wind
xout1nowind = stateout1nowind[:, 0]
zout1nowind = stateout1nowind[:, 1]

m_1st_stage_off = m0 + m_prelaunch + mdotF1 * (delta_t1a + delta_t1b) - mfairing
m_2nd_stage_off = m_1st_stage_off + delta_t2*mdot2  # m0 - m_1st_stage_off + (delta_t2*mdot2)

# Checking the masses after MECO and SECO
print("after 1st stage is off at 156 s: ", m_1st_stage_off)
print("after 2nd stage is off at 529 s: ", m_2nd_stage_off)


def downrange():
    down1 = np.arctan(zout1[-1]/xout1[-1]) * Rplanet
    down1nowind = np.arctan(zout1nowind[-1]/xout1nowind[-1]) * Rplanet
    print(down1/1000, down1nowind/1000)

    return down1, down1nowind


#############################################
############ Plots and grpahs ###############

# print colored trajectory

f_trajectory = open("trajectory.txt", "r")
data_traj = f_trajectory.read().splitlines()
time_traj = []
x_traj = []
z_traj = []


# Altitude
def plotflight():
    fig, (ax1, ax2) = plt.subplots(1, 2)
    ax1.plot(tout,altitude)
    ax1.plot(tout1, altitude1)

    ax1.grid()
    # 2D orbit
    F_stage_counter = 0
    F_sep_stage_counter = 0
    S_stage_counter = 0
    S_sep_stage_counter = 0
    T_stage_counter = 0
    T_sep_stage_counter = 0

    # 1st stage
    for i in tout:
        if i <= delta_t1a+delta_t1b:
            F_stage_counter = F_stage_counter + 1
            F_sep_stage_counter = F_sep_stage_counter + 1
            S_stage_counter = S_stage_counter + 1
            S_sep_stage_counter = S_sep_stage_counter + 1
            T_stage_counter = T_stage_counter + 1
            T_sep_stage_counter = T_sep_stage_counter + 1
        if i > delta_t1a+delta_t1b and i <= t2start:
            F_sep_stage_counter = F_sep_stage_counter + 1
            S_stage_counter = S_stage_counter + 1
            S_sep_stage_counter = S_sep_stage_counter + 1
            T_stage_counter = T_stage_counter + 1
            T_sep_stage_counter = T_sep_stage_counter + 1
        if i > t2start and i <= t2end:
            S_stage_counter = S_stage_counter + 1
            S_sep_stage_counter = S_sep_stage_counter + 1
            T_stage_counter = T_stage_counter + 1
            T_sep_stage_counter = T_sep_stage_counter + 1
        if i > t2end:
            S_sep_stage_counter = S_sep_stage_counter + 1
            T_stage_counter = T_stage_counter + 1
            T_sep_stage_counter = T_sep_stage_counter + 1


    ax2.plot(xout[1:F_stage_counter],zout[1:F_stage_counter], 'g-')
    ax2.plot(xout[F_stage_counter:F_sep_stage_counter],zout[F_stage_counter:F_sep_stage_counter], 'm-')
    ax2.plot(xout[F_sep_stage_counter:S_stage_counter],zout[F_sep_stage_counter:S_stage_counter], 'r-')
    ax2.plot(xout[S_stage_counter:S_sep_stage_counter],zout[S_stage_counter:S_sep_stage_counter], 'm-')
    ax2.plot(xout[S_sep_stage_counter:T_stage_counter],zout[S_sep_stage_counter:T_stage_counter], 'b-')
    ax2.plot(xout[T_stage_counter:T_sep_stage_counter],zout[T_stage_counter:T_sep_stage_counter], 'y-')

    ax2.plot(xout[0],zout[0], 'g*')
    # adding 1st stage falldown:
    ax2.plot(xout1,zout1, 'y-')
    ax2.plot(xout1nowind, zout1nowind, 'y-')
    theta = np.linspace(0, 2*np.pi, 10000)
    xplanet = Rplanet*np.sin(theta)
    zplanet = Rplanet*np.cos(theta)
    xtargetorbit = (Rplanet + 500000) * np.sin(theta)
    ztargetorbit = (Rplanet + 500000) * np.cos(theta)
    ax2.plot(xplanet,zplanet, 'b-')
    ax2.plot(xtargetorbit, ztargetorbit, color='gray')
    ax2.grid()
    ax2.legend()


def plotmap():
    langle = 85 * np.pi / 180
    lpadlon = 27.954896
    lpadlat = 42.085804
    # electronlpadlon = 177.865080
    # electronlpadlat = -39.261626

    down1, down1nowind = downrange()
    distance1 = lpadlon + (down1*np.sin(langle) / Rplanet) * (180 / np.pi) / np.cos(lpadlat * np.pi/180)
    distance1nowind = lpadlon + (down1nowind / Rplanet) * (180 / np.pi) / np.cos(lpadlat * np.pi/180)
    mapbox_access_token = open("mapbox_publictoken.txt").read()

    lat1 = lpadlat + (down1 * np.cos(langle))/ Rplanet * (180 / np.pi)# / np.cos(27.57 * np.pi/180)# / Rplanet * (180 / np.pi) / np.cos(42.085804 * np.pi/180
    lat3 = lpadlat + (down1nowind * np.cos(langle))/ Rplanet * (180 / np.pi)
    lon1 = lpadlon + (down1*np.sin(langle) / Rplanet) * (180 / np.pi) #/ (np.cos(lat1) * np.pi/180)
    lon3 = lpadlon + (down1nowind * np.sin(langle) / Rplanet) * (180 / np.pi)

    fig = go.Figure(go.Scattermapbox(
            lat=[lat3],
            lon=[lon3],
            mode='markers',
            marker=go.scattermapbox.Marker(
                size=16,
                color='rgb(255, 0, 0)'
            ),
            text=['Stage1_nowind'],
        ))
    fig = fig.add_trace(go.Scattermapbox(
        lat=[lpadlat],
        lon=[lpadlon],
        mode='markers',
        marker=go.scattermapbox.Marker(
            size=10,
            color='rgb(0, 0, 255)'
        ),
        text=['LaunchPad'],
    ))


    fig = fig.add_trace(go.Scattermapbox(
        lat=[lat1],
        lon=[lon1],
        mode='markers',
        marker=go.scattermapbox.Marker(
            size=16,
            color='rgb(255, 255, 0)',
            opacity=0.8
        ),
        text=['Stage1_wind'],
    ))

    fig.update_layout(
        hovermode='closest',
        mapbox=dict(
            accesstoken=mapbox_access_token,
            bearing=0,
            center=go.layout.mapbox.Center(
                lat=42,
                lon=23
            ),
            pitch=0,
            zoom=5
        )
    )

    fig.show()




def velocityplot():
    plt.figure()
    plt.plot(tout,velout, label="Ракета")
    #plt.plot(tout1+150, velout1, label="Първа степен")
    plt.xlabel('Време [s]')
    plt.ylabel('Скорост [m/s]')
    plt.grid()
    plt.legend()

# mass
def massplot():
    plt.figure()
    plt.plot(tout,massout)
    plt.xlabel('Време [s]')
    plt.ylabel('Маса [kg]')
    plt.grid()

def heightplot():
    plt.figure()
    plt.plot(tout, altitude, label="Ракета")
    # plt.plot(tout1+150, altitude1, label="Първа степен")
    plt.xlabel('Време [s]')
    plt.ylabel('Височина [m]')
    # plt.legend()
    plt.grid()

def orbitplot():
    # 2D orbit
    plt.figure()
    F_stage_counter = 0
    F_sep_stage_counter = 0
    S_stage_counter = 0
    S_sep_stage_counter = 0
    T_stage_counter = 0
    T_sep_stage_counter = 0

    # 1st stage
    for i in tout:
        if i <= delta_t1a + delta_t1b:
            F_stage_counter = F_stage_counter + 1
            F_sep_stage_counter = F_sep_stage_counter + 1
            S_stage_counter = S_stage_counter + 1
            S_sep_stage_counter = S_sep_stage_counter + 1
            T_stage_counter = T_stage_counter + 1
            T_sep_stage_counter = T_sep_stage_counter + 1
        if i > delta_t1a + delta_t1b and i <= t2start:
            F_sep_stage_counter = F_sep_stage_counter + 1
            S_stage_counter = S_stage_counter + 1
            S_sep_stage_counter = S_sep_stage_counter + 1
            T_stage_counter = T_stage_counter + 1
            T_sep_stage_counter = T_sep_stage_counter + 1
        if i > t2start and i <= t2end:
            S_stage_counter = S_stage_counter + 1
            S_sep_stage_counter = S_sep_stage_counter + 1
            T_stage_counter = T_stage_counter + 1
            T_sep_stage_counter = T_sep_stage_counter + 1
        if i > t2end:
            S_sep_stage_counter = S_sep_stage_counter + 1
            T_stage_counter = T_stage_counter + 1
            T_sep_stage_counter = T_sep_stage_counter + 1


    plt.plot(xout[1:F_stage_counter],zout[1:F_stage_counter], 'g-')
    plt.plot(xout[F_stage_counter:F_sep_stage_counter],zout[F_stage_counter:F_sep_stage_counter], 'k-')
    plt.plot(xout[F_sep_stage_counter:S_stage_counter],zout[F_sep_stage_counter:S_stage_counter], 'r-')
    plt.plot(xout[S_stage_counter:S_sep_stage_counter],zout[S_stage_counter:S_sep_stage_counter], 'm-', label="Орбитален полет")
    plt.plot(xout[S_sep_stage_counter:T_stage_counter],zout[S_sep_stage_counter:T_stage_counter], 'b-')
    plt.plot(xout[T_stage_counter:T_sep_stage_counter],zout[T_stage_counter:T_sep_stage_counter], 'y-')

    plt.plot(xout[0],zout[0], 'g*')
    # adding 1st stage falldown:
    # plt.plot(xout1,zout1, 'y-', label="Първа степен, падане")
    # plt.plot(xout1nowind, zout1nowind, 'y-')
    theta = np.linspace(0, 2*np.pi, 10000)
    xplanet = Rplanet*np.sin(theta)
    zplanet = Rplanet*np.cos(theta)
    xtargetorbit = (Rplanet + 500000) * np.sin(theta)
    ztargetorbit = (Rplanet + 500000) * np.cos(theta)
    plt.plot(xplanet,zplanet, 'b-', label="Земя")
    plt.plot(xtargetorbit, ztargetorbit, color='gray', label="Целева орбита")
    plt.grid()
    plt.legend(loc=3)

plotflight()
velocityplot()
heightplot()
orbitplot()
plotmap()
plt.show()