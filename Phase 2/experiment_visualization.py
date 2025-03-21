# -*- coding: utf-8 -*-
"""
Created on Thu Mar 20 19:19:27 2025

@author: johnt
"""

from scipy.interpolate import interp1d
import numpy as np
import matplotlib.pyplot as plt

"""Define Function for Dictionaries"""

def experiment1():
    
    experiment = {'time_range' : np.array([0,20000]),
                  'initial_conditions' : np.array([0.3025,0]), #[0] rover velocity, [1] rover position
                  'alpha_dist' : np.array([0, 100, 200, 300, 400, 500, 600, \
                                           700, 800, 900, 1000]),
                  'alpha_deg' : np.array([11.509, 2.032, 7.182, 2.478, \
                                        5.511, 10.981, 5.601, -0.184, \
                                        0.714, 4.151, 4.042]),
                  'Crr' : 0.1}
    
    
    # Below are default values for example only:
    end_event = {'max_distance' : 50,
                 'max_time' : 5000,
                 'min_velocity' : 0.01 #incase the rover gets stuck
                 }
    
    return experiment, end_event

"""Initialize Dictionaries"""
experiment, end_event = experiment1()
alpha_dist = experiment["alpha_dist"]
alpha_deg = experiment["alpha_deg"]


alpha_fun = interp1d(alpha_dist, alpha_deg, kind="cubic", fill_value="extrapolate") #Fits the cubic spline")

steps = np.linspace(alpha_dist[0], alpha_dist[-1], 100)
interp_deg = alpha_fun(steps)

plt.plot(steps, interp_deg)
plt.scatter(alpha_dist, alpha_deg)
plt.xlabel("Distance (meters)")
plt.ylabel("Slope Angle (degrees)")
plt.show()


