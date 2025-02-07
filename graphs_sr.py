# -*- coding: utf-8 -*-
"""
Created on Fri Feb  7 14:15:10 2025

@author: johnt
"""

import matplotlib.pyplot as plt
import subfunctions
import numpy as np


rover = subfunctions.rover
get_mass = subfunctions.get_mass
torque = subfunctions.tau_dcmotor
gearbox = subfunctions.get_gear_ratio

gear_ratio = gearbox(rover["wheel_assembly"]["speed_reducer"])

omega_input_shaft = np.linspace(0.01,rover["wheel_assembly"]["motor"]["speed_noload"],num=50)

torque_input_shaft = torque(omega_input_shaft,rover["wheel_assembly"]["motor"])

omega_output_shaft = (1/gear_ratio)*omega_input_shaft

torque_output_shaft = gear_ratio*torque_input_shaft

power_shaft = omega_output_shaft*torque_output_shaft

plt.figure(figsize=(10,10))

"""Subplot for shaft speed vs shaft torque"""
plt.subplot(3,1,1)
plt.plot(torque_output_shaft,omega_output_shaft,marker="o")
plt.xlabel("Motor Shaft Torque (Nm)")
plt.ylabel("Motor Shaft Speed (rad/s)")


"""Subplot for motor power vs shaft torque"""
plt.subplot(3,1,2)
plt.plot(torque_output_shaft,power_shaft,marker="o")
plt.xlabel("Motor Shaft Torque (Nm)")
plt.ylabel("Motor Power (W)")

"""Subplot for motor power vs shaft speed"""
plt.subplot(3,1,3)
plt.plot(omega_output_shaft,power_shaft,marker="o")
plt.xlabel("Motor Shaft Speed (rad/s)")
plt.ylabel("Motor Power (W)")

plt.show()