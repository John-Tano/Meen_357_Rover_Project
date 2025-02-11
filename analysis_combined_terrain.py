"""
Created on monday 2/10/2025 1:33:23 2025

@author: Alier
"""
import numpy as np
import math
import matplotlib.pyplot as plt
import subfunctions
from mpl_toolkits.mplot3d import Axes3D
Crr_array = np.linspace(0.01,0.4,25)
slope_array_deg = np.linspace(-10,35,25)
CRR,Slope =np.meshgrid(Crr_array,slope_array_deg)
VMAX=np.zeros(np.shape(CRR),dtype=float)
def bisection(fun, lb, ub, err_max, iter_max,Crr,planet,rover,angle):

    import numpy as np
    if not callable(fun):
        raise Exception("First argument must be a callable function.")
    if not (isinstance(lb, (int, float)) and isinstance(ub, (int, float))):
        raise Exception("Lower and upper bounds must be scalar values.")
    if not (lb<ub):
        raise Exception("Lower has to be less.")
    if not (isinstance(err_max, (int, float)) and err_max > 0):
        raise Exception("Maximum error must be a positive scalar.")
    if not (isinstance(iter_max, int) and iter_max > 0):
        raise Exception("Maximum iterations must be a positive integer.")
    

    if np.isnan(fun(lb,angle,rover,planet,Crr)) or np.isnan(fun(ub,angle,rover,planet,Crr)) or np.isinf(fun(ub,angle,rover,planet,Crr)) or np.isinf(fun(ub,angle,rover,planet,Crr)):
        return None, None, 0, -2
    if fun(lb,angle,rover,planet,Crr) * fun(ub,angle,rover,planet,Crr)> 0:
        return None, None, 0, -1
    
    numIter = 0
    err = float('inf')
    root = lb  
    f_lb=fun(lb,angle,rover,planet,Crr)
    while numIter < iter_max and err > err_max:
        prev_root = root
        root = (lb + ub) / 2  # Midpoint
        f_root = fun(root,angle,rover,planet,Crr)

        
        if np.isnan(f_root) or np.isinf(f_root):
            return None, None, numIter, -2
        
        # Compute relative error (except on first iteration)
        if numIter > 0:
            err = abs((root - prev_root) / root)*100
        
        # Bisection logic
        if f_root == 0 or err < err_max:
            return root, err, numIter, 1  # Root found with sufficient accuracy
        elif f_lb * f_root < 0:
            ub = root  # Root is in lower subinterval
        else:
            lb = root  # Root is in upper subinterval
            f_lb = f_root
        
        numIter += 1
    
    exitFlag = 1 if err < err_max else 0
    return root, err, numIter, exitFlag
v_maxes =[]
err_max=1e-5
N=np.shape(CRR)[0]
iter_max=15
fun=subfunctions.F_net
lb=0
ub=subfunctions.rover['wheel_assembly']['motor']['speed_noload']
planet=subfunctions.planet
rover=subfunctions.rover
for i in range(N):
    for j in range(N):
        Crrpart=float(CRR[i,j])
        slope2 =float(Slope[i,j])
        maxo=(bisection(fun,0,rover['wheel_assembly']['motor']['speed_noload'],err_max,iter_max,Crrpart,planet,rover,slope2)[0])
        if maxo!=None:
            VMAX[i][j]=float(maxo)/subfunctions.get_gear_ratio(rover['wheel_assembly']['speed_reducer'])* rover['wheel_assembly']['wheel']['radius']
        else:
            VMAX[i][j]=None
figures=plt.figure()
ax=Axes3D(figures,elev=20,azim=20)
ax.set_xlabel("Rolling resitance Coefficent")
ax.set_ylabel("Degree of Angle")
ax.set_zlabel("velocity max")
ax.plot_surface(CRR,Slope,VMAX)
ax.set_title("velocity max vs Degree of Angle and Rolling resitance Coefficent ")
figures.add_axes(ax)
plt.show()