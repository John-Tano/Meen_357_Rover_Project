import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import root_scalar
from subfunctions import *
w = np.linspace(0,3.8,25)
slope = np.linspace(-10,35,25)
Crr = 0.2
gear = get_gear_ratio(rover["wheel_assembly"]["speed_reducer"])
r = rover["wheel_assembly"]["wheel"]["radius"]
wn = rover["wheel_assembly"]["motor"]["speed_noload"]

wmax = np.zeros(len(slope), dtype = float)
v_max = np.zeros(len(slope), dtype = float)

for i in range(len(wmax)):
    function = lambda w: F_net(w, float(slope[i]), rover, planet, Crr)
    solution = root_scalar(function, method='bisect', bracket=[0, wn])
    wmax[i] = solution.root  
    v_max[i] = (wmax[i]) / (gear*r)

plt.xlabel("Terrain Angle (deg)")
plt.ylabel("Max Rover Speed (m/s)")
plt.title("Max Rover Speed vs Terrain Angle")
plt.plot(slope, v_max, label = 'Max Speed')
plt.grid(True)
plt.legend()
plt.show()
