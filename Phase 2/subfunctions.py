
import math
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





"""Initializing Dictionaries"""

rover, planet = define_rover_1()
experiment, end_event = experiment1()
alpha_dist = experiment["alpha_dist"]
alpha_deg = experiment["alpha_deg"]




"""Phase 1"""

def get_mass(rover):
    """
    Inputs:  rover:  dict      Data structure containing rover parameters
    
    Outputs:     m:  scalar    Rover mass [kg].
    """
    
    # Check that the input is a dict
    if type(rover) != dict:
        raise Exception('Input must be a dict')
    
    # add up mass of chassis, power subsystem, science payload, 
    # and components from all six wheel assemblies
    m = rover['chassis']['mass'] \
        + rover['power_subsys']['mass'] \
        + rover['science_payload']['mass'] \
        + 6*rover['wheel_assembly']['motor']['mass'] \
        + 6*rover['wheel_assembly']['speed_reducer']['mass'] \
        + 6*rover['wheel_assembly']['wheel']['mass'] \
    
    return m


def get_gear_ratio(speed_reducer):
    """
    Inputs:  speed_reducer:  dict      Data dictionary specifying speed
                                        reducer parameters
    Outputs:            Ng:  scalar    Speed ratio from input pinion shaft
                                        to output gear shaft. Unitless.
    """
    
    # Check that the input is a dict
    if type(speed_reducer) != dict:
        raise Exception('Input must be a dict')
    
    # Check 'type' field (not case sensitive)
    if speed_reducer['type'].lower() != 'reverted':
        raise Exception('The speed reducer type is not recognized.')
    
    # Main code
    d1 = speed_reducer['diam_pinion']
    d2 = speed_reducer['diam_gear']
    
    Ng = (d2/d1)**2
    
    return Ng


def tau_dcmotor(omega, motor):
    """
    Inputs:  omega:  numpy array      Motor shaft speed [rad/s]
             motor:  dict             Data dictionary specifying motor parameters
    Outputs:   tau:  numpy array      Torque at motor shaft [Nm].  Return argument
                                      is same size as first input argument.
    """
    
    # Check that the first input is a scalar or a vector
    if (type(omega) != int) and (type(omega) != float) and (not isinstance(omega, np.ndarray)):
        raise Exception('First input must be a scalar or a vector. If input is a vector, it should be defined as a numpy array.')
    elif not isinstance(omega, np.ndarray):
        omega = np.array([omega],dtype=float) # make the scalar a numpy array
    elif len(np.shape(omega)) != 1:
        raise Exception('First input must be a scalar or a vector. Matrices are not allowed.')

    # Check that the second input is a dict
    if type(motor) != dict:
        raise Exception('Second input must be a dict')
        
    # Main code
    tau_s    = motor['torque_stall']
    tau_nl   = motor['torque_noload']
    omega_nl = motor['speed_noload']
    
    # initialize
    tau = np.zeros(len(omega),dtype = float)
    for ii in range(len(omega)):
        if omega[ii] >= 0 and omega[ii] <= omega_nl:
            tau[ii] = tau_s - (tau_s-tau_nl)/omega_nl *omega[ii]
        elif omega[ii] < 0:
            tau[ii] = tau_s
        elif omega[ii] > omega_nl:
            tau[ii] = 0
        
    return tau
    
    


def F_rolling(omega, terrain_angle, rover, planet, Crr):
    """
    Inputs:           omega:  numpy array     Motor shaft speed [rad/s]
              terrain_angle:  numpy array     Array of terrain angles [deg]
                      rover:  dict            Data structure specifying rover 
                                              parameters
                    planet:  dict            Data dictionary specifying planetary 
                                              parameters
                        Crr:  scalar          Value of rolling resistance coefficient
                                              [-]
    
    Outputs:           Frr:  numpy array     Array of forces [N]
    """
    
    # Check that the first input is a scalar or a vector
    if (type(omega) != int) and (type(omega) != float) and (not isinstance(omega, np.ndarray)):
        raise Exception('First input must be a scalar or a vector. If input is a vector, it should be defined as a numpy array.')
    elif not isinstance(omega, np.ndarray):
        omega = np.array([omega],dtype=float) # make the scalar a numpy array
    elif len(np.shape(omega)) != 1:
        raise Exception('First input must be a scalar or a vector. Matrices are not allowed.')
        
    # Check that the second input is a scalar or a vector
    if (type(terrain_angle) != int) and (type(terrain_angle) != float) and (not isinstance(terrain_angle, np.ndarray)):
        raise Exception('Second input must be a scalar or a vector. If input is a vector, it should be defined as a numpy array.')
    elif not isinstance(terrain_angle, np.ndarray):
        terrain_angle = np.array([terrain_angle],dtype=float) # make the scalar a numpy array
    elif len(np.shape(terrain_angle)) != 1:
        raise Exception('Second input must be a scalar or a vector. Matrices are not allowed.')
        
    # Check that the first two inputs are of the same size
    if len(omega) != len(terrain_angle):
        raise Exception('First two inputs must be the same size')
    
    # Check that values of the second input are within the feasible range  
    if max([abs(x) for x in terrain_angle]) > 75:    
        raise Exception('All elements of the second input must be between -75 degrees and +75 degrees')
        
    # Check that the third input is a dict
    if type(rover) != dict:
        raise Exception('Third input must be a dict')
        
    # Check that the fourth input is a dict
    if type(planet) != dict:
        raise Exception('Fourth input must be a dict')
        
    # Check that the fifth input is a scalar and positive
    if (type(Crr) != int) and (type(Crr) != float):
        raise Exception('Fifth input must be a scalar')
    if Crr <= 0:
        raise Exception('Fifth input must be a positive number')
        
    # Main Code
    m = get_mass(rover)
    g = planet['g']
    r = rover['wheel_assembly']['wheel']['radius']
    Ng = get_gear_ratio(rover['wheel_assembly']['speed_reducer'])
    
    v_rover = r*omega/Ng
    
    Fn = np.array([m*g*math.cos(math.radians(x)) for x in terrain_angle],dtype=float) # normal force
    Frr_simple = -Crr*Fn # simple rolling resistance
    
    Frr = np.array([math.erf(40*v_rover[ii]) * Frr_simple[ii] for ii in range(len(v_rover))], dtype = float)
    
    return Frr


def F_gravity(terrain_angle, rover, planet):
    """
    Inputs:  terrain_angle:  numpy array   Array of terrain angles [deg]
                     rover:  dict          Data structure specifying rover 
                                            parameters
                    planet:  dict          Data dictionary specifying planetary 
                                            parameters
    
    Outputs:           Fgt:  numpy array   Array of forces [N]
    """
    
    # Check that the first input is a scalar or a vector
    if (type(terrain_angle) != int) and (type(terrain_angle) != float) and (not isinstance(terrain_angle, np.ndarray)):
        raise Exception('First input must be a scalar or a vector. If input is a vector, it should be defined as a numpy array.')
    elif not isinstance(terrain_angle, np.ndarray):
        terrain_angle = np.array([terrain_angle],dtype=float) # make the scalar a numpy array
    elif len(np.shape(terrain_angle)) != 1:
        raise Exception('First input must be a scalar or a vector. Matrices are not allowed.')
        
    # Check that values of the first input are within the feasible range  
    if max([abs(x) for x in terrain_angle]) > 75:    
        raise Exception('All elements of the first input must be between -75 degrees and +75 degrees')

    # Check that the second input is a dict
    if type(rover) != dict:
        raise Exception('Second input must be a dict')
    
    # Check that the third input is a dict
    if type(planet) != dict:
        raise Exception('Third input must be a dict')
        
    # Main Code
    m = get_mass(rover)
    g = planet['g']
    
    Fgt = np.array([-m*g*math.sin(math.radians(x)) for x in terrain_angle], dtype = float)
        
    return Fgt


def F_drive(omega, rover):
    """
    Inputs:  omega:  numpy array   Array of motor shaft speeds [rad/s]
             rover:  dict          Data dictionary specifying rover parameters
    
    Outputs:    Fd:  numpy array   Array of drive forces [N]
    """
    
    # Check that the first input is a scalar or a vector
    if (type(omega) != int) and (type(omega) != float) and (not isinstance(omega, np.ndarray)):
        raise Exception('First input must be a scalar or a vector. If input is a vector, it should be defined as a numpy array.')
    elif not isinstance(omega, np.ndarray):
        omega = np.array([omega],dtype=float) # make the scalar a numpy array
    elif len(np.shape(omega)) != 1:
        raise Exception('First input must be a scalar or a vector. Matrices are not allowed.')

    # Check that the second input is a dict
    if type(rover) != dict:
        raise Exception('Second input must be a dict')
    
    # Main code
    Ng = get_gear_ratio(rover['wheel_assembly']['speed_reducer'])
    
    tau = tau_dcmotor(omega, rover['wheel_assembly']['motor'])
    tau_out = tau*Ng
    
    r = rover['wheel_assembly']['wheel']['radius']
    
    # Drive force for one wheel
    Fd_wheel = tau_out/r 
    
    # Drive force for all six wheels
    Fd = 6*Fd_wheel
    
    return Fd


def F_net(omega, terrain_angle, rover, planet, Crr):
    """
    Inputs:           omega:  list     Motor shaft speed [rad/s]
              terrain_angle:  list     Array of terrain angles [deg]
                      rover:  dict     Data structure specifying rover 
                                      parameters
                     planet:  dict     Data dictionary specifying planetary 
                                      parameters
                        Crr:  scalar   Value of rolling resistance coefficient
                                      [-]
    
    Outputs:           Fnet:  list     Array of forces [N]
    """
    
    # Check that the first input is a scalar or a vector
    if (type(omega) != int) and (type(omega) != float) and (not isinstance(omega, np.ndarray)):
    # if (not isinstance(omega, np.ndarray)):
        raise Exception('First input must be a scalar or a vector. If input is a vector, it should be defined as a numpy array.')
    elif not isinstance(omega, np.ndarray):
        omega = np.array([omega],dtype=float) # make the scalar a numpy array
    elif len(np.shape(omega)) != 1:
        raise Exception('First input must be a scalar or a vector. Matrices are not allowed.')
        
    # Check that the second input is a scalar or a vector
    if (type(terrain_angle) != int) and (type(terrain_angle) != float) and (not isinstance(terrain_angle, np.ndarray)):
        raise Exception('Second input must be a scalar or a vector. If input is a vector, it should be defined as a numpy array.')
    elif not isinstance(terrain_angle, np.ndarray):
        terrain_angle = np.array([terrain_angle],dtype=float) # make the scalar a numpy array
    elif len(np.shape(terrain_angle)) != 1:
        raise Exception('Second input must be a scalar or a vector. Matrices are not allowed.')
        
    # Check that the first two inputs are of the same size
    if len(omega) != len(terrain_angle):
        raise Exception('First two inputs must be the same size')
    
    # Check that values of the second input are within the feasible range  
    if max([abs(x) for x in terrain_angle]) > 75:    
        raise Exception('All elements of the second input must be between -75 degrees and +75 degrees')
        
    # Check that the third input is a dict
    if type(rover) != dict:
        raise Exception('Third input must be a dict')
        
    # Check that the fourth input is a dict
    if type(planet) != dict:
        raise Exception('Fourth input must be a dict')
        
    # Check that the fifth input is a scalar and positive
    if (type(Crr) != int) and (type(Crr) != float):
        raise Exception('Fifth input must be a scalar')
    if Crr <= 0:
        raise Exception('Fifth input must be a positive number')
    
    # Main Code
    Fd = F_drive(omega, rover)
    Frr = F_rolling(omega, terrain_angle, rover, planet, Crr)
    Fg = F_gravity(terrain_angle, rover, planet)
    
    Fnet = Fd + Frr + Fg # signs are handled in individual functions
    
    return Fnet

omega = 1 # rad/s (motor shaft speed)
angle = 5 # degrees (terrain angle)
Crr = 0.1

Fd = F_drive(omega, rover)
Frr = F_rolling(omega, angle, rover, planet, Crr)
Fg = F_gravity(angle, rover, planet)
Fnet = F_net(omega, angle, rover, planet, Crr)

"""PHASE 2 OF THE PROJECT"""

def motorW(v, rover):
    """
    Computes the shaft speed of motor in rad/s
    
    Inputs:
        v: 1D numpy array or scalar, rover translation velocity (m/s)
        rover: dict, data structure containing rover definition

    
    Returns:
        w: 1D numpy aray or scalar, motor shaft speed (rad/s)
    """
    #Check arguments for correct data type
    if not (isinstance(v,(int,float)) or (isinstance(v,np.ndarray) and v.ndim==1)):
        raise Exception("The translational velocity argument must be a 1-D array or scalar")
    if not isinstance(rover, dict):
        raise Exception("The rover argument must be a dictionary")
        
    #Initialize motor speed based on data type of translational velocity
    if isinstance(v,(int,float)):
        w = 0
    if isinstance(v,dict):
        w = np.zeros(v)
        
    ng = get_gear_ratio(rover["wheel_assembly"]["speed_reducer"])
    final_ng = ng/(rover["wheel_assembly"]["wheel"]["radius"])
    
    w = v*final_ng
    return w

def rover_dynamics(t, y, rover, planet, experiment):
    """
    Computes the derivative of the state vector for the rover.
    
    Inputs:
        t: scalar, time sample [s]
        y: 1D numpy array, state vector [velocity (m/s), position (m)]
        rover: dict, data structure containing rover definition
        planet: dict, data structure containing planet definition
        experiment: dict, data structure containing experiment definition
    
    Returns:
        dydt: 1D numpy array, derivative of state vector [acceleration (m/s^2), velocity (m/s)]
    """
    # Input validation
    if not (isinstance(t, (float, int)) or isinstance(t, np.ndarray)):
        raise Exception("time argument must be a scalar or numpy array")
    if not (isinstance(y, np.ndarray) and y.ndim == 1 and len(y) == 2):
        raise Exception("state vector argument must be a 2 element numpy array")
    if not (isinstance(rover, dict) and isinstance(planet, dict) and isinstance(experiment, dict)):
        raise Exception("rover, planet, and experiment arguments must be dictionaries")
        
    # Interpolate terrain angle at the current rover position
    alpha_fun = interp1d(alpha_dist, alpha_deg, kind="cubic", fill_value="extrapolate")
    terrain_angle = alpha_fun(y[1])  # Terrain angle at current position [deg]
    
    # Current velocity and mass
    v_0 = y[0]  # Rover velocity [m/s]
    m = get_mass(rover)  # Rover mass [kg]
    
    # Calculate motor angular velocity
    w_m = motorW(float(v_0), rover)  # Motor angular velocity [rad/s]
    
    # Calculate net force
    crr = experiment["Crr"]  # Rolling resistance coefficient
    F = F_net(w_m, float(terrain_angle), rover, planet, crr)  # Net force [N]
    
    # Calculate acceleration
    acceleration = F / m  # Acceleration [m/s^2]
    
    # Return derivative of state vector
    dydt = np.array([acceleration.item(), v_0])  # [acceleration (m/s^2), velocity (m/s)]
    return dydt

def mechpower(v,rover):
    """
    Computes instantaneous power output of a sing;e motor based on translation velocity
    
    Inputs:
        v: 1D numpy array or scalar, rover translation velocity [m/s]
        rover: dict, data structure containing rover definition
        
    Returns:
        P: 1D numpy array or scalar, instantaneous power output of a single motor based on velocity
    """
    
    if not (isinstance(v,(int,float)) or (isinstance(v,np.ndarray) and v.ndim==1)):
        raise Exception("The translational velocity argument must be a 1-D array or scalar")
    if not isinstance(rover, dict):
        raise Exception("The rover argument must be a dictionary")
        
    # if isinstance(v,(int,float)):
    #     P = 0
    # if (isinstance(v,np.ndarray) and v.ndim==1):
    #     P = np.zeros(v)
    
    motor = rover["wheel_assembly"]["motor"]
    
    omega_motor = motorW(v,rover)
    tau_motor = tau_dcmotor(omega_motor,motor)
    
    P = tau_motor * omega_motor

    return P

def battenergy(t, v, rover):
    """
    Computes total energy consumed from the battery pack by the rover
    
    Inputs:
        t: 1D numpy array, Array of time samples from a rover simulation [s]
        v: 1D numpy array or scalar, rover translation velocity [m/s]
        rover: dict, data structure containing rover definition
        
    Returns:
        E: scalar, Total electrical energy consumed from the rover battery pack over the simulation [J]
    """
    if not (isinstance(t,np.ndarray) and t.ndim==1):
        raise Exception("Time argument must be a 1-D numpy array")
    if not (isinstance(v,np.ndarray) and v.ndim==1):
        raise Exception("The translational velocity argument must be a 1-D array or scalar")
    if not isinstance(rover, dict):
        raise Exception("The rover argument must be a dictionary")
    if not (len(t)==len(v)):
        raise Exception("Time and velocity arguments must be the same length numpy arrays")
        
    effcy = rover["wheel_assembly"]["motor"]["effcy"]
    effcy_tau = rover["wheel_assembly"]["motor"]["effcy_tau"]
    omega_motor = motorW(v,rover)

    tau_motor = tau_dcmotor(omega_motor,rover)
    
    effcy_fun = interp1d(effcy_tau, effcy, kind="cubic")
    interp_effcy = effcy_fun(tau_motor) #Motor efficiency as a function of torque in decimal form
    p_instant = (mechpower(v,rover)) / interp_effcy #Array of instantaneous power consumption accounting for efficiency at each time/ velocity
    E = np.trapz(p_instant,t)
    

    return E
















