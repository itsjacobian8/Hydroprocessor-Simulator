#!/usr/bin/env python3
from time import ctime, time
from termcolor import colored

class GetInputs:
    """
    collection of miscellaneous methods for obtaining user inputs
    
    class attributes
    ----------------
    please note this class is not instantiated with any inputs.
    
    the following attributes are assigned after the class has been instantiated
        inputs: dictionary for storing user inputs
        system: string specifying the system to be simulated (bubble column vs ebullated bed)
        recycleInfo: string specifying if the system has no, constant or variable internal recycle
    
    class methods
    -------------
    getSystem: None
        allows user to specify the system to be simulated
    
    getRecycleInfo: None
        allows user to specify if the system has no, constant or variable internal recycle
    
    getInletFlowRates: None
        allows the user to specify inlet gas and liquid flow rates
    
    getGasLiquidFlowRates: None
        allows user to specify gas-liquid properties and parameters
    
    getSolidProperties: None
        allows the user to specify solid properties and parameters
    
    getGridParameters: None
        allows the user to specify distributor grid parameters
    
    getColumnParameters: None
        allows the user to specify column geometrical parameters and operating conditions
    
    getConvergenceCriteria: None
        allows user to specify convergence criteria
    
    for more information, consult docstrings for individual methods
    """
    def __init__(self):
        self.inputs = {}
        self.system = None
        self.recycleInfo = None
        
    def getSystem(self):
        """
        allows user to specify which system to simulate
        """
        system = str(input(f"[{ctime(time())}] choose one of: 'bubble column', 'ebullated bed': "))
        
        if(system.startswith("bubble") and system.endswith("column")):
            print(colored(f"[{ctime(time())}] You've selected 'bubble column.'\n", "green", "on_black"))
        elif(system.startswith("ebullated") and system.endswith("bed")):
            print(colored(f"[{ctime(time())}] You've selected 'ebullated bed.'\n", "green", "on_black"))
        else:
            raise Exception("Please enter one of the options exactly as provided on the prompt!")
        
        self.system = system
        
    def getRecycleInfo(self):
        """
        allows user to specify internal recycle for bubble column or ebullated bed
        """
        if(self.system.startswith("bubble") and self.system.endswith("column")):
            print(colored(f"[{ctime(time())}] Which of the following describes your internal recycle situation?", "green", "on_black"))
            print(colored(f"[{ctime(time())}] 1. no recycle", "white", "on_black"))
            print(colored(f"[{ctime(time())}] 2. constant recycle", "cyan", "on_black"))
            recycle = str(input(f"[{ctime(time())}] Please enter one string exactly as shown above: "))
            if(recycle.startswith("no") and recycle.endswith("recycle")):
                print(colored(f"[{ctime(time())}] You've selected 'no recycle'\n", "green", "on_black"))
            elif(recycle.startswith("constant") and recycle.endswith("recycle")):
                print(colored(f"[{ctime(time())}] You've selected 'constant recycle'\n", "green", "on_black"))
            else:
                raise Exception("Enter one option exactly as shown on prompt!")
        elif(self.system.startswith("ebullated") and self.system.endswith("bed")):
            print(colored(f"[{ctime(time())}] Which of the following describes your internal recycle situation?", "green", "on_black"))
            print(colored(f"[{ctime(time())}] 1. no recycle", "white", "on_black"))
            print(colored(f"[{ctime(time())}] 2. variable recycle", "cyan", "on_black"))
            recycle = str(input(f"[{ctime(time())}] Please enter one string exactly as shown above: "))
            if(recycle.startswith("no") and recycle.endswith("recycle")):
                print(colored(f"[{ctime(time())}] You've selected 'no recycle' \n", "green", "on_black"))
            elif(recycle.startswith("variable") and recycle.endswith("recycle")):
                print(colored(f"[{ctime(time())}] You've selected 'variable recycle' \n", "green", "on_black"))
            else:
                raise Exception("Enter one option exactly as shown on prompt!")
        
        self.recycleInfo = recycle
        
    def getInletFlowRates(self):
        """
        obtains inlet flow rates from user
        """
        self.inputs['qgtMin'] = float(input(f"[{ctime(time())}] Enter minimum inlet gas flow rate [m^3/s]: "))
        print(colored(f"[{ctime(time())}] Minimum inlet gas flow = {self.inputs['qgtMin']}\n","green", "on_black"))
        
        self.inputs['qgtMax'] = float(input(f"[{ctime(time())}] Enter maximum inlet gas flow rate [m^3/s]: "))
        print(colored(f"[{ctime(time())}] Maximum inlet gas flow rate = {self.inputs['qgtMax']}\n", "green", "on_black"))
        
        self.inputs['qgtIncrement'] = float(input(f"[{ctime(time())}] Enter inlet gas flow rate increment [m^3/s]: "))
        print(colored(f"[{ctime(time())}] Inlet gas flow rate increment = {self.inputs['qgtIncrement']}\n", "green", "on_black"))
        
        self.inputs['qlvr'] = float(input(f"[{ctime(time())}] Enter inlet liquid flow rate [m^3/s]: "))
        print(colored(f"[{ctime(time())}] Inlet liquid flow rate = {self.inputs['qlvr']}\n", "green", "on_black"))
    
    def getGasLiquidProperties(self):
        """
        obtains gas and liquid parameters from user
        """
        self.inputs['rhog'] = float(input(f"[{ctime(time())}] Enter gas density [kg/m^3]: "))
        print(colored(f"[{ctime(time())}] Gas density = {self.inputs['rhog']}\n", "green", "on_black"))
        
        self.inputs['rhol'] = float(input(f"[{ctime(time())}] Enter liquid density [kg/m^3]: "))
        print(colored(f"[{ctime(time())}] Liquid density = {self.inputs['rhol']}\n", "green", "on_black"))
        
        self.inputs['mul'] = float(input(f"[{ctime(time())}] Enter liquid viscosity [Pa.s]: "))
        print(colored(f"[{ctime(time())}] Liquid viscosity = {self.inputs['mul']}\n", "green", "on_black"))
        
        self.inputs['sigma'] = float(input(f"[{ctime(time())}] Enter surface tension [N/m]: "))
        print(colored(f"[{ctime(time())}] Surface tension = {self.inputs['sigma']}\n", "green", "on_black"))
        
        self.inputs['classWidth'] = float(input(f"[{ctime(time())}] Enter bubble class width [mm]: "))
        print(colored(f"[{ctime(time())}] Bubble class width = {self.inputs['classWidth']}\n", "green", "on_black"))
        
        self.inputs['alpha'] = float(input(f"[{ctime(time())}] Enter confidence interval for bubble size distribution: "))
        print(colored(f"[{ctime(time())}] Confidence interval = {self.inputs['alpha']}\n", "green", "on_black"))
        
    def getSolidProperties(self):
        """ 
        obtains solid phase parameters from user 
        """
        self.inputs['rhop'] = float(input(f"[{ctime(time())}] Enter particle density [kg/m^3]: "))
        print(colored(f"[{ctime(time())}] Particle density = {self.inputs['rhop']}\n", "green", "on_black"))
        
        self.inputs['Vp'] = float(input(f"[{ctime(time())}] Enter particle volume [m^3]: "))
        print(colored(f"[{ctime(time())}] Particle volume = {self.inputs['Vp']}\n", "green", "on_black"))
        
        self.inputs['Ap'] = float(input(f"[{ctime(time())}] Enter particle surface area [m^2]: "))
        print(colored(f"[{ctime(time())}] Particle surface area = {self.inputs['Ap']}\n", "green", "on_black"))
        
        self.inputs['mcat'] = float(input(f"[{ctime(time())}] Enter mass of catalyst inventory [kg]: "))
        print(colored(f"[{ctime(time())}] Mass of catalyst inventory {self.inputs['mcat']}\n", "green", "on_black"))
        
        self.inputs['hset'] = float(input(f"[{ctime(time())}] Enter desired bed/column height [m]: "))
        print(colored(f"[{ctime(time())}] Desired bed/column height = {self.inputs['hset']}\n", "green", "on_black"))
        
    def getGridParameters(self):
        """
        obtains distributor grid parameters
        """
        self.inputs['dor'] = float(input(f"[{ctime(time())}] Enter outlet orifice diameter [m]: "))
        print(colored(f"[{ctime(time())}] Outlet orifice diameter = {self.inputs['dor']}\n", "green", "on_black"))
        
        self.inputs['nor'] = float(input(f"[{ctime(time())}] Enter number of outlet orifices per riser: "))
        print(colored(f"[{ctime(time())}] Number of outlet orifices = {self.inputs['nor']}\n", "green", "on_black"))
        
        self.inputs['nr'] = float(input(f"[{ctime(time())}] Enter number of risers on the distributor: "))
        print(colored(f"[{ctime(time())}] Number of risers = {self.inputs['nr']}\n", "green", "on_black"))
        
    def getColumnParameters(self):
        """
        obtains column parameters
        """
        self.inputs['P'] = float(input(f"[{ctime(time())}] Enter operating pressure [MPa]: "))
        print(colored(f"[{ctime(time())}] Operating pressure = {self.inputs['P']}\n", 'green', "on_black"))
        
        self.inputs['Dc'] = float(input(f"[{ctime(time())}] Enter column diameter [m]: "))
        print(colored(f"[{ctime(time())}] Column diameter = {self.inputs['Dc']}\n", "green", "on_black"))
        
        self.inputs['drec'] = float(input(f"[{ctime(time())}] Enter recycle line diameter [m]: "))
        print(colored(f"[{ctime(time())}] Recycle line diameter = {self.inputs['drec']}\n", "green", "on_black"))
        
        self.inputs['vsep'] = float(input(f"[{ctime(time())}] Enter gas-liquid separator volume [m^3]: "))
        print(colored(f"[{ctime(time())}] Gas-liquid separator volume = {self.inputs['vsep']}\n", "green", "on_black"))
        
    def getConvergenceCriteria(self):
        """
        obtains convergence criteria
        """
        if(self.recycleInfo.startswith("constant") and self.recycleInfo.endswith("recycle")):
            self.inputs['R'] = float(input(f"[{ctime(time())}] Enter constant liquid recycle ratio: "))
            print(colored(f"[{ctime(time())}] Constant liquid recycle ratio = {self.inputs['R']}\n", "green", "on_black"))
            
        self.inputs['Rstep'] = float(input(f"[{ctime(time())}] Enter liquid recycle ratio step size: "))
        print(colored(f"[{ctime(time())}] Liquid recycle ratio step size = {self.inputs['Rstep']}\n", "green", "on_black"))
        
        self.inputs['absTol'] = float(input(f"[{ctime(time())}] Enter absolute tolerance: "))
        print(colored(f"[{ctime(time())}] Absolute tolerance = {ctime(time())}\n", "green", "on_black"))
        
        self.inputs['relTol'] = float(input(f"[{ctime(time())}] Enter relative tolerance: "))
        print(colored(f"[{ctime(time())}] Relative tolerance = {self.inputs['relTol']} \n", "green", "on_black"))
        
        self.inputs['maxIter'] = int(input(f"[{ctime(time())}] Enter maximum number of iterations: "))
        print(colored(f"[{ctime(time())}] Maximum number of iterations = {self.inputs['maxIter']}\n", "green", "on_black"))
