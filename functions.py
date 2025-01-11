import numpy as np
import scipy as sc
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import matplotlib.animation as anim

year=31536000 ##nb de seconde en 1 année

rR= 140000 ##140 000 km 
Mp=1
Msat=5.68e26
G=6.6743e-11*Msat*year**2/1e9 ##conversion en masse de saturne (pour avoir directement q et non la masse) et année
Rp=60268 #60 268 km
D=3.7e-4
TR=2*np.pi/np.sqrt(G*Mp/rR**3) ##conversion en masse de saturne et année
Tdisk=3.1*10**5 
K=39*2.3e-4*Rp**5*(G*Mp)**(1/2)/2
F=23.5*D**3/TR


def eq_DELTA(delta,t,M,Mp):
    delta_p=2**5/3**3*M/Mp*D/TR*delta**(-3)
    return delta_p

def eq_R(delta,t,M,Mp):
    r=rR*(delta+1)
    r_p=3*Rp**5*2.3e-4*(G*Mp)**(1/2)*M/Mp*r**(-11/2)
    return r_p


class planet:

    def __init__(self, dt):
        self.r = rR
        self.delta = 1e-8
        self.M = F*dt
        self.regim = "continue"
        self.age = 0
        self.transition = 0

      
    
    def iteration(self, system, index_planet, nb_accretion, dt):
        
        self.age += dt
        def continue_regim(self, system, nb_accretion):
            self.M += F*dt
            self.delta = odeint(eq_DELTA,self.delta,np.linspace(self.age,self.age+dt,100),args=(self.M, Mp))[-1][0] 
            + (odeint(eq_R,self.r,np.linspace(self.age,self.age+dt,100),args=(self.M,Mp))[-1][0] -rR)/rR
            self.r = rR*(self.delta+1)

            if self.delta >= 8.4*D:

                self.regim = "discret"
                self.transition = self.age

            return nb_accretion

        def discret_regim(self, system, index_planet, nb_accretion):

            i = index_planet

            self.delta = odeint(eq_DELTA,self.delta,np.linspace(self.age,self.age+dt,100),args=(self.M, Mp))[-1][0] 
            + (odeint(eq_R,self.r,np.linspace(self.age,self.age+dt,100),args=(self.M,Mp))[-1][0] -rR)/rR

            self.r = rR*(self.delta+1)

            if len(system)>1 and i != len(system)-1 :
                delta_accretion = 2*((system[i].M+system[i+1].M)/(3*Mp))**(1/3)
                if abs(system[i].delta-system[i+1].delta) <= delta_accretion:
                    system[i].M = system[i].M + system[i+1].M 
                    system.pop(i+1)
                    nb_accretion= nb_accretion+ 1

            return nb_accretion


        
        if self.regim == "continue":

            nb_accretion = continue_regim(self, system, nb_accretion)

        if self.regim == "discret": 

            nb_accretion = discret_regim(self, system, index_planet, nb_accretion)
        
        return nb_accretion
    


        






