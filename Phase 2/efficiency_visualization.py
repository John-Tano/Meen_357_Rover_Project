# -*- coding: utf-8 -*-
"""
Created on Sat Mar 22 03:56:25 2025

@author: johnt
"""

import numpy as np
from scipy.interpolate import interp1d
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
 
"""Dictionaries Section """ 

def define_rover_1():
    """ Initialize Rover dict for testing """
    wheel = {'radius':0.30,
             'mass':1}
    speed_reducer = {'type':'reverted',
                     'diam_pinion':0.04,
                     'diam_gear':0.07,
                     'mass':1.5}
    
    motor = {'torque_stall':170,
             'torque_noload':0,
             'speed_noload':3.80,
             'mass':5.0,
             'effcy_tau':np.array([0, 10, 20, 40, 70, 165]),
             'effcy':np.array([0, 0.55, 0.75, 0.71, 0.50, 0.05])
             }
    # telemetry = {
    #     'Time': ,
    #     'completion_time': ,
    #     'velocity': ,
    #     'position': ,
    #     'distance_traveled: ,
    #     'max_velocity': ,
    #     'average_velocity': ,
    #     'power': ,
    #     'batery_energy': ,
    #     'energy_per_distance': 
    #     }
    
    
        
    chassis = {'mass':659}
    science_payload = {'mass':75}
    power_subsys = {'mass':90}

    wheel_assembly = {'wheel':wheel,
                      'speed_reducer':speed_reducer,
                      'motor':motor}
    
    rover = {'wheel_assembly':wheel_assembly,
             'chassis':chassis,
             'science_payload':science_payload,
             'power_subsys':power_subsys
             # 'telemetry':telemetry
             }
    
    planet = {'g':3.72}
    
    # return everything we need
    return rover, planet

rover, planet = define_rover_1()

# torque_stall = rover["wheel_assembly"]["motor"]["torque_stall"]
# torque_noload = rover["wheel_assembly"]["motor"]["torque_noload"]
effcy = rover["wheel_assembly"]["motor"]["effcy"]
effcy_tau = rover["wheel_assembly"]["motor"]["effcy_tau"]

steps = np.linspace(effcy_tau[0], effcy_tau[-1], 100) #100 evenly spaced torque values from 0 to 165 N*M

effcy_fun = interp1d(effcy_tau, effcy, kind="cubic")
interp_effcy = 100 * effcy_fun(steps) #Motor efficiency in terms of percent

plt.figure()
plt.plot(steps,interp_effcy)
plt.scatter(effcy_tau, effcy*100, color = "red", marker = "*")
plt.xlabel("Motor Torque (N*M)") 
plt.ylabel("Motor Efficiency (%)")
plt.show()






























