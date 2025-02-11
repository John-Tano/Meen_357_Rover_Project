# -*- coding: utf-8 -*-
"""
Created on Wed Feb  5 09:48:34 2025

@author: johnt
"""
import matplotlib.pyplot as plt
import subfunctions
import numpy as np


rover = subfunctions.rover
get_mass = subfunctions.get_mass
torque = subfunctions.tau_dcmotor

omega_shaft = np.linspace(0.01,rover["wheel_assembly"]["motor"]["speed_noload"],num=50)

torque_shaft = torque(omega_shaft,rover["wheel_assembly"]["motor"])

power_shaft = omega_shaft*torque_shaft

plt.figure(figsize=(10,10))

"""Subplot for shaft speed vs shaft torque"""
plt.subplot(3,1,1)
plt.plot(torque_shaft,omega_shaft,marker="o")
plt.xlabel("Motor Shaft Torque (Nm)")
plt.ylabel("Motor Shaft Speed (rad/s)")


"""Subplot for motor power vs shaft torque"""
plt.subplot(3,1,2)
plt.plot(torque_shaft,power_shaft,marker="o")
plt.xlabel("Motor Shaft Torque (Nm)")
plt.ylabel("Motor Power (W)")

"""Subplot for motor power vs shaft speed"""
plt.subplot(3,1,3)
plt.plot(omega_shaft,power_shaft,marker="o")
plt.xlabel("Motor Shaft Speed (rad/s)")
plt.ylabel("Motor Power (W)")
plt.tight_layout()
plt.show()
