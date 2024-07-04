import numpy as np
# constants
G = 6.6742*10**-11 # gravitational constant (SI unit)

# Planet Earth
Rplanet = 6371000 # m
Mplanet = 5.972e24 # kg

# Wind
alpha = 90 * np.pi / 180
windfactor = 1.0

# Parameters of rocket
D = 1.2 # meters Diameter of rocket
S = np.pi*(D/2)**2  # cross sectional area of rocket
L1 = D * 12.1  # windward area 1st stage
L2 = D * 2.4  # windward area 2nd stage
CD = 0.4
cdwind = 0.82  # for a long cylinder

#initial conditions for single stage rocket
x0 = Rplanet
z0 = 0.0
velx0 = 0.0
velz0 = 0.0
period = 1200.0
period1 = 800

# Mass
m1all = 10200 + 50 # kg
m1prop = 9250 # kg

m2all = 2400  # kg
m2prop = 2150  # kg

mfairing = 50  # kg

m3all = 150  # kg

m0 = m1all + m2all + m3all # kg The mass of the whole rocket

#thrust
thrust1 = 192000.0 # N
thrust2 = 18980.0  # N

# Specific Impulse
Isp1 = 292.0
Isp2 = 350.0

mdotF1 = - thrust1 / (Isp1 * 9.81)
mdot2 = - thrust2 / (Isp2 * 9.81)


# Duration Times in seconds
# MECO 149s, 1st stage separation at 153s, 2nd stage ignition at 156s, 2nd stage cut off at 529s
t0 = 3.0
t_vert = 12.0 # working 12.0
t_pitchover = 60.0 # working 60.0 #85.0
delta_t1a = 140.0
delta_t1b = 9.0
delta_t1off = 4.0
t1end = delta_t1a + delta_t1b + delta_t1off
delta_t02 = 3.0
t2start = t1end + delta_t02
t2man = 190.0
delta_t2 = 373.0
delta_t20ff = 1.0
t2end = t2start + delta_t2

psi1 = 0.0
dpsi1 = 22.0

psi2 = 0.0
dpsi2 = 0.0

psi3 = 0.0
dpsi3 = 0.0

psi4 = 80.0
dpsi4 = 0.0

psi5 = 0.0
dpsi5 = 0.0



thrustfactor1 = 1
thrustfactor2 = 0.8
thrustfactor3 = 1
thrustfactor4 = 1
thrustfactor5 = 1
thrustfactor6 = 1

m_prelaunch = t0 * mdotF1
m_launch = m0 + m_prelaunch
mdry1 = m1all + m_prelaunch + mdotF1 * (delta_t1a + delta_t1b) - mfairing




