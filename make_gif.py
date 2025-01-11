import numpy as np
import scipy as sc
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import matplotlib.animation as anim
import functions as fc
import importlib
import os
from PIL import Image


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
dt = 1e-3
t=np.arange(0,1,dt)


planets = [fc.planet(dt)]

for i in range(len(t)):
    nb_accretion = 0


    for j in range(len(planets)):
        j = j - nb_accretion ###probleme avec ca
        nb_accretion = planets[j].iteration(planets, j,nb_accretion, dt)
        
    if planets[-1].delta >= 8.4*D:
        planets.append(fc.planet(dt))
    plt.figure()
    for k in range(len(planets)):

        plt.scatter(planets[k].r, 0, s= 3e4*(planets[k].M)**(1/3))

    plt.xlim([140000,150000])
    plt.ylim([-0.02,0.02])
    plt.axis("off")
    plt.savefig(f"./figures/time_{i}.png")
    plt.close()


images = [Image.open(f"./figures/time_{i}.png") for i in range(len(t))]

images[0].save('accretion.gif', save_all=True, append_images=images[1:], duration=50, loop=0)