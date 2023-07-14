sim_type = "planar" # planar or spherical

"""Parameters"""
T = 0.592 # kcal/mol
kc = 20. # bending rigidity: kcal/mol
chi = 0.4 # coverage
c0 = 340. # Spontaneous curvature preference of lipid: microns^-1
r = 8.442806489 # vesicle radius in microns
D0 = 8.0 # diffusion constant: microns^2/s
Ap = 6.34e-7 # Area per lipid: microns^2
eta = 0.128 # water viscotiy at 25C: kcal s/mol microns^3
mode_max = 1.0 #30 # max mode, l for SH and q for planar

"""Time parameters"""
rate = 1./100. # rate of saving positions
Nsteps = 7000 # number of steps

