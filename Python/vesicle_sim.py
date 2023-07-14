#!/opt/local/bin/python
import numpy as np
import sys
from sys import argv
from scipy.optimize import curve_fit
from scipy.special import lpmn, factorial

"""
FOR USE WITH PYTHON 3
TO RUN: python vesicle_sim.py parameter_file uq_output_file pq_output_file 
"""

filename = argv[1]
exec(open(filename).read())

uq_filename = argv[2]
f1 = open(uq_filename, 'w')
pq_filename = argv[3]
f2 = open(pq_filename, 'w')

"""Calculated Parameters"""
A = 4*np.pi*r**2
if chi == 0:
    c0 = 0.0
    Ap == 0:
    rho0 = 0.
else:
    rho0 = 2.*chi/Ap

if sim_type == "planar":
    q_max = mode_max
    P_max = int((q_max*r)-2)
    P = np.arange(0,P_max+1)
    Q = 0
    q = (P+2)/r
    Lambda = 1/(4 * eta * q) # Hydrodynamic tensor for membrane equation of motion                        
    if chi == 0:
        D = 0.0
    else:
        D = rho0*D0*q**2 # diffusion constant used in equation of motion
    taum = 4. * eta / (kc * q**3) # relaxation time of membrane                                                     
    taup = 1. / (q**2 * D0) # relaxation time of particle distribution                                                

    """define empty arrays for complex uq and pq"""
    uq = np.zeros((q.shape[0]),dtype=complex)
    pq = np.zeros((q.shape[0]),dtype=complex)

if sim_type == "spherical":
    l_max = mode_max
    l = np.arange(2,l_max+1)
    m = np.repeat(np.arange(0,l_max+1),2)[1:]
    Lambda = (1/(eta*r**3)) * (l*(l+1.)/((2.*l+1.)*(2.*l**2.+2.*l-1.))) # Hydrodynamic tensor for membrane
    if chi == 0:
        D = 0.0
    else:
        D = (2*l*(l+1)*chi*(1-chi)*D0)/(Ap*r**4) # diffusion constant used in equation of motion
    taum = (1/(kc*Lambda*l*(l+1)*(l+2)*(l-1))) # relaxation time of membrane                                         
    taup = r**2 / (l*(l+1) * D0) # relaxation time for the particle distribution                                      

    """define empty arrays for ulm and rholm"""
    ulm = np.zeros((len(l), len(m)))
    rholm = np.zeros((len(l), len(m)))

    """constants for projection into equitorial plane"""
    plm = lpmn(m[-1],l_max,np.cos(np.pi/2))[0].T
    plm = np.repeat(plm,2,axis=1)[2:,1:]
    Nlm = np.sqrt((2*l[:,None]+1)*factorial(l[:,None]-m[None,:])/(4*np.pi*factorial(l[:,None]+m[None,:])))

"""Time parameters"""
dt = (taum[-1])/100 # time step                                                                                        
tmax = Nsteps*rate / dt # max amount of time to run

"""Main loop"""
step = 0

while step < tmax:
    if sim_type == "planar":
        randm = np.sqrt(T*dt*Lambda*A) * (np.random.normal(size=q.shape[0])+1j*np.random.normal(size=q.shape[0]))
        randp = np.sqrt(q**2*D*dt*A) * (np.random.normal(size=q.shape[0])+1j*np.random.normal(size=q.shape[0]))
        
        fuq = - kc * q**4 * uq + Ap * c0 * (kc/2.) * q**2 * pq
        if chi == 0:
            fpq = 0.0
        else:
            fpq = - Ap*pq*T/(2.*chi*(1-chi)) + Ap * c0 * (kc/2.) * q**2 * uq
            
        uq += dt * Lambda * fuq + randm
        pq += dt * D * fpq / T + randp
        if step % int(rate/dt) == 0:
            print(step, end = " ", file=f1)
            for i in range(P.shape[0]):
                print(0, uq[i].real, 0, uq[i].imag, end=" ", file=f1)
            print("\n", file=f1)
            print(step, end = " ", file=f2)
            for i in range(P.shape[0]):
                print(0, pq[i].real, 0, pq[i].imag, end=" ", file=f2)
            print("\n", file=f2)

    if sim_type == "spherical":
        randm = np.sqrt(2*T*dt*Lambda[:,None]) * np.random.normal(size=(len(l),len(m)))
        randp = np.sqrt(2*D[:,None]*dt) * np.random.normal(size=(len(l),len(m)))
            
        fulm = - (kc * (l[:,None]+2) * (l[:,None]-1) * (l[:,None] *(l[:,None]+1)) * ulm + kc*c0*Ap*r*(-2+l[:,None]+l[:,None]**2)*rholm/2.)
        if chi == 0:
            frholm = 0.0
        else:
            frholm = - ((T*Ap*r**2)/(2*chi*(1 - chi)) * rholm + kc*c0*Ap*r*(-2+l[:,None]+l[:,None]**2)*ulm/2.)
            
        ulm += dt * Lambda[:,None] * fulm + randm
        rholm += dt * (D[:,None]/T) * frholm + randp
        vq = np.sum(ulm*Nlm*plm,axis=0)
        rhoq = np.sum(rholm*Nlm*plm, axis=0)
        if step % int(rate/dt) == 0:
            print(step, end=" ", file=f1)
            for i in range(3,len(m)-1,2):
                print(0, vq[i], 0, vq[i+1], end=" ", file=f1)
            print("\n", file=f1)
            print(step, end=" ", file=f2)
            for i in range(3,len(m)-1,2):
                print(0, rhoq[i], 0, rhoq[i+1], end=" ", file=f2)
            print("\n", file=f2)

    step += 1


