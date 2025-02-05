# -*- coding: utf-8 -*-
"""
Created on Fri Jan 24 11:31:32 2025

@author: johnt
"""
import numpy as np
import scipy
from scipy.special import erf

rover = {
    
    "wheel_assembly" : {
        "wheel" : {
            "radius":0.30,#wheel radius (m) 
            "mass":1#mass of one drive wheel (kg)
          },
        "speed_reducer" : {
            "type":"reverted", #the type of gearbox
            "diam_pinion":0.04,#speed reducer pinion diameter (m)
            "diam_gear":0.07,#speed reducer gear diameter (m)
            "mass":1.5, #speed reducer mass (kg)
            },
        "motor" : {
            "torque_stall":170,#motor stall torque in (Nm)
            "torque_noload":0,#motor no-load torque (Nm)
            "speed_noload":3.8,#motor no-load speed (rad/s)
            "mass":5#motor mass (kg)
            }},
    
    
    "chassis" : {
        "mass":659#chassis mass (kg)
        },
    
    
    "science_payload" : {
        "mass":75#science payload mass (kg)
        },
    
    
    "power_subsys" : {
        "mass":90 #power subsystem mass (kg)
        },
    }

planet = {
    "g":3.72#acceleration due to gravity (m/s^2)
    }

def get_mass(rover):
    if not isinstance(rover,dict):
        raise Exception("input must be a dictionary")
    total_mass=0
    
    for key, component in rover.items():
        if isinstance(component, dict):
            for subkey, subcomponent in component.items():
                if isinstance(subcomponent,dict):
                    for sskey,sscomponent in subcomponent.items():
                        if sskey == "mass":
                            total_mass += sscomponent*6
                if subkey == "mass":
                    total_mass += subcomponent #currently not accessing masses in wheel assembly
    return total_mass

def get_gear_ratio(speed_reducer):
    if not isinstance(speed_reducer,dict):
        raise ValueError("input must be a dictionary")
    Ng=0
    d_pinion=0
    d_gear=0
    
    for key, component in speed_reducer.items():
        if key=="type":
            if not component.lower() == "reverted":
                raise ValueError("The type of speed reducer must be reverted")
        if key == "diam_pinion":
            d_pinion=component
        if key == "diam_gear":
            d_gear=component
    if d_gear==0 or d_pinion==0:
        raise ValueError("gear and pinion diameters need to be nonzero")
    Ng = (d_gear/d_pinion)**2
    return Ng

def tau_dcmotor(omega, motor):
    if not (isinstance(omega, (int,float)) or (isinstance(omega, np.ndarray) and omega.ndim==1)):
        raise ValueError("angular velocity must be either a 1D numpy array or scalar")
    if not isinstance(motor, dict):
        raise ValueError("the motor argument must be passed as a dicitonary")
    t_S=0
    t_Nl=0
    w_Nl=0
    
    for key, component in motor.items():
        if key == "torque_stall":
            t_S=component
        if key == "torque_noload":
            t_Nl=component
        if key == "speed_noload":
            w_Nl=component
    
    if isinstance(omega,(int,float)):
        if omega>w_Nl:
            return 0
        elif omega<0:
            return t_S
        else:
            return t_S-((t_S-t_Nl)/w_Nl)*omega
    if isinstance(omega, np.ndarray):
        tau=np.zeros_like(omega)
        for i,o in enumerate(omega):
            if o>w_Nl:
                tau[i]=0
            elif o<0:
                tau[i]=t_S
            else:
                tau[i] = t_S-((t_S-t_Nl)/w_Nl)*o
    return tau

def F_drive(omega, rover):
    if not (isinstance(omega,np.ndarray) or (isinstance(omega, (int,float)) and abs(omega) <=75)):
        raise ValueError("angular velocity must be either an array or scalar")
    if not isinstance(rover, dict):
        raise ValueError("the motor argument must be passed as a dicitonary")
    
    if isinstance(omega,(int,float)):
        F_onewheel = tau_dcmotor(omega,rover["wheel_assembly"]["motor"])*get_gear_ratio(rover["wheel_assembly"]["speed_reducer"])/(rover["wheel_assembly"]["wheel"]["radius"])
        
    if isinstance(omega, np.ndarray):
        F_onewheel = tau_dcmotor(omega,rover["wheel_assembly"]["motor"])*get_gear_ratio(rover["wheel_assembly"]["speed_reducer"])/(rover["wheel_assembly"]["wheel"]["radius"])
    Fd=F_onewheel*6
    
    return Fd

def F_gravity(terrain_angle, rover, planet):
    if not ((isinstance(terrain_angle,np.ndarray) and np.all(abs(terrain_angle)<=75)) or isinstance(terrain_angle,(float,int))):
        raise ValueError("terrain _angle argument must be a numpy array and all angles be less than 75 degrees")
    if not isinstance(rover, dict):
        raise ValueError("the motor argument must be passed as a dicitonary")
    if not isinstance(planet, dict):
        raise ValueError("the plant argument must be passed as a dicitonary")
   
    rad_angle = np.radians(terrain_angle)
    m_rover = get_mass(rover)
    g = planet["g"]
    
    Fgt=-1*m_rover*g*np.sin(rad_angle)

    return Fgt

def F_rolling(omega, terrain_angle,rover,planet,Crr):
    if not (isinstance(omega,np.ndarray) or isinstance(omega, (int,float))):
        raise ValueError("angular velocity must be either an array or scalar")
    if not ((isinstance(terrain_angle,np.ndarray) and np.all(abs(terrain_angle)<=75)) or isinstance(terrain_angle,(float,int))):
        raise ValueError("terrain _angle argument must be a numpy array and all angles be less than 75 degrees")
        
    if not (((isinstance(omega,np.ndarray) and isinstance(terrain_angle,np.ndarray)) and omega.size == terrain_angle.size) or (isinstance(omega,(float,int)) and isinstance(terrain_angle,(float,int)))):
        raise ValueError("angular velocity and terrain_angle arguments must both be scalars or same length numpy arrays")
    if not isinstance(rover, dict):
        raise ValueError("the rover argument must be passed as a dicitonary")
    if not isinstance(planet, dict):
        raise ValueError("the planet argument must be passed as a dicitonary")
    if not (isinstance(Crr,(int,float)) and Crr>0):
        raise ValueError("rolling resistance coefficient argument must be a positive scalar")
     
    rad_angle = np.radians(terrain_angle)
    m_rover = get_mass(rover)
    g = planet["g"]
    r_W = rover["wheel_assembly"]["wheel"]["radius"]
    Ng = get_gear_ratio(rover["wheel_assembly"]["speed_reducer"])
    omega_out = omega/Ng
    
    v_rover = r_W*omega_out
    F_Norm = m_rover*g*np.cos(rad_angle)
    Frr_simple = Crr*F_Norm
    Frr_corrected = -1*erf(40*v_rover)*Frr_simple
    
    return Frr_corrected

def F_net(omega, terrain_angle, rover, planet, Crr):
    if not (isinstance(omega,np.ndarray) or isinstance(omega, (int,float))):
        raise ValueError("angular velocity must be either an array or scalar")
    if not ((isinstance(terrain_angle,np.ndarray) and np.all(abs(terrain_angle)<=75)) or isinstance(terrain_angle,(float,int))):
        raise ValueError("terrain _angle argument must be a numpy array and all angles be less than 75 degrees")
        
    if not (((isinstance(omega,np.ndarray) and isinstance(terrain_angle,np.ndarray)) and omega.size == terrain_angle.size) or (isinstance(omega,(float,int)) and isinstance(terrain_angle,(float,int)))):
        raise ValueError("angular velocity and terrain_angle arguments must both be scalars or same length numpy arrays")
    if not isinstance(rover, dict):
        raise ValueError("the rover argument must be passed as a dicitonary")
    if not isinstance(planet, dict):
        raise ValueError("the planet argument must be passed as a dicitonary")
    if not (isinstance(Crr,(int,float)) and Crr>0):
        raise ValueError("rolling resistance coefficient argument must be a positive scalar")
      
    Fd = F_drive(omega, rover)
    Fgt = F_gravity(terrain_angle,rover, planet)
    Frr = F_rolling(omega,terrain_angle,rover,planet,Crr)
    
    Fn=Fd+Fgt+Frr
    
    return Fn


