import numpy as np
# constants
G = 6.6742*10**-11 # gravitational constant (SI unit)

# Planet Earth
Rplanet = 6371000 # m
Mplanet = 5.972e24 # kg

# Wind
# wind0 = 10  # m/s
alpha = 90 * np.pi / 180
windfactor = 1.0

# Parameters of rocket
D = 10.0 # meters Diameter of rocket
S = np.pi*(D/2)**2  # cross sectional area of rocket
L1 = D * 42  # windward area 1st stage
L2 = D * 24  # windward area 2nd stage
CD = 0.4
cdwind = 0.82  # for a long cylinder

#initial conditions for single stage rocket
x0 = Rplanet
z0 = 0.0
velx0 = 0.0
velz0 = 0.0
period = 12000
period1 = 360
period2 = 870

# Mass
# m0 = 2822000 # kg The mass of the whole rocket

m1all = 2214000 # kg
m1prop = 2077000 # kg

m2all = 470000 # kg
m2prop = 456100 # kg

m3all = 120000 # kg
m3prop = 107800 # kg

m0 = m1all + m2all + m3all # kg The mass of the whole rocket

#thrust
thrustF1 = 6672*10**3 # N single engine - 5xthrustF1
thrustJ2 = 1023*10**3 # N
thrust2 = 5*thrustJ2 # N
thrust3 = 1*thrustJ2 # N

# Specific Impulse
Isp1 = 260.0  # seconds
Isp2 =421.0  # seconds
Isp3 = 421.0  # seconds

mdotF1 = - thrustF1 / (Isp1 * 9.81)
mdot2 = - thrust2 / (Isp2 * 9.81)
mdot3 = - thrust3 / (Isp3 * 9.81)

# Duration Times in seconds
t0 = 6.0
t_vert = 12.0
t_pitchover = 100.0
delta_t1a = 135.0
delta_t1b = 15.0
delta_t1off = 1.0
t1end = delta_t1a + delta_t1b + delta_t1off
delta_t02 = 1.0
t2start = t1end + delta_t02
t2man = 300.0
delta_t2 = 365.0
delta_t20ff = 3.0
t2end = t2start + delta_t2
t3start = t2end + delta_t20ff
delta_t3 = 142.0
t3end = t3start + delta_t3

psi1 = 0.0
dpsi1 = 30.0

psi2 = 0.0
dpsi2 = 0.0

psi3 = 0.0
dpsi3 = 0.0

psi4 = 24.6
dpsi4 = 0.0

psi5 = 0.0
dpsi5 = 0.0

psi6 = 0.0
dpsi6 = 0.0


thrustfactor1 = 0.95
thrustfactor2 = 0.95
thrustfactor3 = 0.98
thrustfactor4 = 1
thrustfactor5 = 0.98
thrustfactor6 = 1

m_prelaunch = t0 * 5 * mdotF1
m_launch = m0 + m_prelaunch





