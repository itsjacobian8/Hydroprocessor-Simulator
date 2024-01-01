#!/usr/bin/env python3
from time import time, ctime
from termcolor import colored

def getSystem():
    """
    allows user to specify which system they wish to simulate: bubble column vs ebullated bed
    
    return:
        system: switch for bubble column or ebullated bed (boolean)
    """
    system = str(input(f"[{ctime(time())}] Choose one of: 'bubble column', 'ebullated bed': "))
    bubbleColumn = (system.startswith("bubble") and system.endswith("column"))
    ebullatedBed = (system.startswith("ebullated") and system.endswith("bed"))
    
    if(bubbleColumn):
        print(colored(f"[{ctime(time())}] You've selected 'bubble column'\n", "green"))
    elif(ebullatedBed):
        print(colored(f"[{ctime(time())}] You've selected 'ebullated bed'\n", "green"))
    else:
        raise Exception("Please enter one of the options exactly as provided on the prompt!")
        
    return system
    
def getRecycleInfo():
    """
    allows user to specify internal recycle for bubble column or ebullated bed
    
    return:
        output: switch for internal recycle (boolean)
    """
    recycleInfo = str(input(f"[{ctime(time())}] Does your system have an internal recycle? (yes/no): "))
    recycle = (recycleInfo == "yes")
    noRecycle = (recycleInfo == "no")
    
    if(recycle):
        output = True
        print(colored(f"You've chosen a system with an internal recycle.\n","green"))
    elif(noRecycle):
        output = False
        print(colored(f"You've chosen a system without an internal recycle.\n","green"))
    else:
        raise Exception("You need to choose one option exactly as shown on the prompt!")
        
    return output
    
def getBubbleColumnParameters(recycleInfo):
    """
    obtains bubble column parameters from user
    
    args:
        recycleInfo: switch for internal recycle (boolean)
    
    return:
        inputs: user supplied input parameters for bubble column (dictionary)
    """
    inputs = {}

    inputs["qgMin"] = float(input(f"[{ctime(time())}] Enter minimum inlet gas flow rate [m^3/s]: "))
    print(colored(f"[{ctime(time())}] Minimum inlet gas flow rate = {inputs['qgMin']} \n","green"))

    inputs["qgMax"] = float(input(f"[{ctime(time())}] Enter maximum inlet gas flow rate [m^3/s]: "))
    print(colored(f"[{ctime(time())}] Maximum inlet gas flow rate = {inputs['qgMax']} \n","green"))

    inputs['qgIncrement'] = float(input(f"[{ctime(time())}] Enter inlet gas flow rate increment [m^3/]: "))
    print(colored(f"[{ctime(time())}] Inlet gas flow rate increment = {inputs['qgIncrement']} \n","green"))
    
    inputs["qlIn"] = float(input(f"[{ctime(time())}] Enter inlet liquid flow rate [m^3/s]: "))
    print(colored(f"[{ctime(time())}] Inlet liquid flow rate = {inputs['qlIn']} \n","green"))
    
    inputs["rhog"] = float(input(f"[{ctime(time())}] Enter gas density [kg/m^3]: "))
    print(colored(f"[{ctime(time())}] Gas density = {inputs['rhog']} \n","green"))
    
    inputs["rhol"] = float(input(f"[{ctime(time())}] Enter liquid density [kg/m^3]: "))
    print(colored(f"[{ctime(time())}] Liquid Density = {inputs['rhol']} \n","green"))
    
    inputs["sigma"] = float(input(f"[{ctime(time())}] Enter gas-liquid interfacial tension [N/m]: "))
    print(colored(f"[{ctime(time())}] Gas-liquid interfacial tension = {inputs['sigma']} \n","green"))
    
    inputs["mul"] = float(input(f"[{ctime(time())}] Enter liquid viscosity [pa.s]: "))
    print(colored(f"[{ctime(time())}] Liquid viscosity = {inputs['mul']} \n","green"))
    
    inputs["P"] = float(input(f"[{ctime(time())}] Enter operating pressure [MPa]: "))
    print(colored(f"[{ctime(time())}] Operating pressure = {inputs['P']} \n","green"))
    
    inputs["Dc"] = float(input(f"[{ctime(time())}] Enter column diameter [m]: "))
    print(colored(f"[{ctime(time())}] Column diameter = {inputs['Dc']} \n","green"))
    
    inputs["dor"] = float(input(f"[{ctime(time())}] Enter outlet orifice diameter [m]: "))
    print(colored(f"[{ctime(time())}] Outlet orifice diameter = {inputs['dor']} \n","green"))
    
    inputs["nor"] = float(input(f"[{ctime(time())}] Enter number of outlet orifices per riser: "))
    print(colored(f"[{ctime(time())}] Number of outlet orifices = {inputs['nor']} \n","green"))
    
    inputs["nr"] = float(input(f"[{ctime(time())}] Enter number of risers on the distributor grid:"))
    print(colored(f"[{ctime(time())}] Number of risers = {inputs['nr']} \n","green"))
    
    inputs["absTol"] = float(input(f"[{ctime(time())}] Enter absolute tolerance: "))
    print(colored(f"[{ctime(time())}] Absolute tolerance = {inputs['absTol']} \n","green"))
    
    inputs["relTol"] = float(input(f"[{ctime(time())}] Enter relative tolerance: "))
    print(colored(f"[{ctime(time())}] Relative tolerance = {inputs['relTol']} \n","green"))
    
    inputs["classWidth"] = float(input(f"[{ctime(time())}] Enter bubble class width [mm]: "))
    print(colored(f"[{ctime(time())}] Bubble class width = {inputs['classWidth']} \n","green"))
    
    inputs["alpha"] = float(input(f"[{ctime(time())}] Enter confidence interval for bubble size distribution: "))
    print(colored(f"[{ctime(time())}] Confidence interval = {inputs['alpha']} \n","green"))
    
    if(recycleInfo):
        inputs["R"] = float(input(f"[{ctime(time())}] Enter liquid recycle ratio :"))
        print(colored(f"[{ctime(time())}] Liquid recycle ratio = {inputs['R']}\n", "green"))
        
        inputs["vsep"] = float(input(f"[{ctime(time())}] Enter separator volume [m^3]:"))
        print(colored(f"[{ctime(time())}] Separator volume = {inputs['vsep']}\n", "green"))
        
        inputs["drec"] = float(input(f"[{ctime(time())}] Enter recycle line diameter [m]:"))
        print(colored(f"[{ctime(time())}] Recycle line diameter = {inputs['drec']} \n", "green"))
        
    return inputs
    
def getEbullatedBedParameters(recycleInfo):
    """
    obtains ebullated bed input parameters from user
    
    args:
        recycleInfo: switch for internal recycle (boolean)
    
    return:
        inputs: user supplied parameters for ebullated bed (dictionary)
    """
    inputs = {}
    inputs["qgtMin"] = float(input(f"[{ctime(time())}] Enter minimum inlet gas flow rate [m^3/s]: "))
    print(colored(f"[{ctime(time())}] Minimum inlet gas flow rate = {inputs['qgtMin']} \n","green"))

    inputs["qgtMax"] = float(input(f"[{ctime(time())}] Enter maximum inlet gas flow rate [m^3/s]: "))
    print(colored(f"[{ctime(time())}] Maximum inlet gas flow rate = {inputs['qgtMax']} \n","green"))

    inputs['qgtIncrement'] = float(input(f"[{ctime(time())}] Enter inlet gas flow rate increment [m^3/]: "))
    print(colored(f"[{ctime(time())}] Inlet gas flow rate increment = {inputs['qgtIncrement']} \n","green"))
    
    inputs["qlvr"] = float(input(f"[{ctime(time())}] Enter inlet liquid flow rate [m^3/s]: "))
    print(colored(f"[{ctime(time())}] Inlet liquid flow rate = {inputs['qlvr']}\n", "green"))
  
    inputs["rhog"] = float(input(f"[{ctime(time())}] Enter gas density [kg/m^3]: "))
    print(colored(f"[{ctime(time())}] Gas density = {inputs['rhog']}\n", "green"))
     
    inputs["rhol"] = float(input(f"[{ctime(time())}] Enter liquid density [kg/m^3]: "))
    print(colored(f"[{ctime(time())}] Liquid density = {inputs['rhol']}\n", "green"))
    
    inputs["mul"] = float(input(f"[{ctime(time())}] Enter liquid viscosity [Pa.s]: "))
    print(colored(f"[{ctime(time())}] Liquid viscosity = {inputs['mul']}\n", "green"))
    
    inputs["mcat"] = float(input(f"[{ctime(time())}] Enter mass of catalyst inventory [kg]: "))
    print(colored(f"[{ctime(time())}] Mass of catalyst inventory = {inputs['mcat']}\n", "green"))
    
    inputs["rhop"] = float(input(f"[{ctime(time())}] Enter particle density [kg/m^3]: "))
    print(colored(f"[{ctime(time())}] Particle density = {inputs['rhop']}\n", "green"))
    
    inputs["Vp"] = float(input(f"[{ctime(time())}] Enter particle volume [m^3]: "))
    print(colored(f"[{ctime(time())}] Particle volume = {inputs['Vp']}\n", "green"))
    
    inputs["Ap"] = float(input(f"[{ctime(time())}] Enter particle surface area [m^2]: "))
    print(colored(f"[{ctime(time())}] Particle surface area = {inputs['Ap']}\n", "green"))
    
    inputs["sigma"] = float(input(f"[{ctime(time())}] Enter gas-liquid interfacial tension [N/m]: "))
    print(colored(f"[{ctime(time())}] Gas-liquid interfacial tension = {inputs['sigma']}\n", "green"))
    
    inputs["P"] = float(input(f"[{ctime(time())}] Enter operating pressure [MPa]: "))
    print(colored(f"[{ctime(time())}] Operating pressure = {inputs['P']}\n", "green"))
    
    inputs["Dc"] = float(input(f"[{ctime(time())}] Enter column diameter [m]: "))
    print(colored(f"[{ctime(time())}] Column Diameter = {inputs['Dc']}\n", "green"))
    
    if(recycleInfo):
        inputs["drec"] = float(input(f"[{ctime(time())}] Enter recycle line diameter [m]: "))
        print(colored(f"[{ctime(time())}] Recycle line diameter = {inputs['drec']}\n", "green"))
    
        inputs["vsep"] = float(input(f"[{ctime(time())}] Enter separator volume [m^3]: "))
        print(colored(f"[{ctime(time())}] Separator volume = {inputs['vsep']}\n", "green"))
        
        inputs["Rstep"] = float(input(f"[{ctime(time())}] Enter liquid recycle ratio step size: "))
        print(colored(f"[{ctime(time())}] Liquid recycle ratio step size = {inputs['Rstep']}\n", "green"))
        
        inputs["maxIter"] = int(input(f"[{ctime(time())}] Enter maximum number of iterations: "))
        print(colored(f"[{ctime(time())}] Maximum number of iterations = {inputs['maxIter']}\n", "green"))
        
    inputs["hset"] = float(input(f"[{ctime(time())}] Enter desired bed height or maximum column height [m]: "))
    print(colored(f"[{ctime(time())}] Desired bed height = {inputs['hset']}\n", "green"))
    
    inputs["dor"] =  float(input(f"[{ctime(time())}] Enter outlet orifice diameter [m]: "))
    print(colored(f"[{ctime(time())}] Outlet orifice diameter = {inputs['dor']}\n", "green"))
    
    inputs["nor"] =  float(input(f"[{ctime(time())}] Enter number of outlet orifices per riser: "))
    print(colored(f"[{ctime(time())}] Number of outlet orifices per riser = {inputs['nor']}\n", "green"))
    
    inputs["nr"] = float(input(f"[{ctime(time())}] Enter number of risers on the distributor grid: "))
    print(colored(f"[{ctime(time())}] Number of risers = {inputs['nr']}\n", "green"))
    
    inputs["classWidth"] = float(input(f"[{ctime(time())}] Enter bubble class width [mm]: "))
    print(colored(f"[{ctime(time())}] Bubble class width = {inputs['classWidth']}\n", "green"))
    
    inputs["alpha"] = float(input(f"[{ctime(time())}] Enter confidence interval for bubble size distribution: "))
    print(colored(f"[{ctime(time())}] Confidence interval = {inputs['alpha']}\n", "green"))
    
    inputs["absTol"] = float(input(f"[{ctime(time())}] Enter absolute tolerance: "))
    print(colored(f"[{ctime(time())}] Absolute tolerance = {inputs['absTol']}\n", "green"))
    
    inputs["relTol"] = float(input(f"[{ctime(time())}] Enter relative tolerance: "))
    print(colored(f"[{ctime(time())}] Relative tolerance = {inputs['relTol']}\n", "green"))
    
    return inputs