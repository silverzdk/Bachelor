import numpy as np 
import scipy as sp
import matplotlib.pyplot as plt

from pyfluids import Fluid, FluidsList, Input

class Pipe:
    def __init__(self,length,outer_diameter,wall_thickness,angle):
        self.length = length
        self.outer_diameter = outer_diameter
        self.wall_thickness = wall_thickness
        self.angle          = angle
        
        self.inner_diameter  = self.outer_diameter - 2 * self.wall_thickness
        self.inner_perimeter = np.pi * self.inner_diameter
        self.outer_perimeter = np.pi * self.outer_diameter

class Node:
    def __init__(self,pipe,x_pos,T_hot_init,T_cold_init,P_hot_init,P_cold_init,m_dot_hot,m_dot_cold,is_boundary=False):
        self.pipe           = pipe
        self.is_boundary    = is_boundary
        self.x_pos          = x_pos

        #Hot flyid
        self.T_hot          = T_hot_init
        self.P_hot          = P_hot_init
        self.m_dot_hot      = m_dot_hot
        
        #cold fluid
        self.T_cold          = T_cold_init
        self.P_cold          = P_cold_init
        self.m_dot_cold      = m_dot_cold

        # Calculate Fluid properties
        
        # Hot fluid
        self.Fluid_hot          = Fluid(FluidsList.Water)
        self.Fluid_hot.update(Input.temperature(self.T_hot), Input.pressure(self.P_hot))
        
        # Hot fluid assuming 100% liquid phase
        self.Fluid_hot_liquid       =Fluid(FluidsList.Water)
        self.Fluid_hot_liquid.update(Input.temperature(self.T_hot),Input.quality(0))

        self.__update_fluid_properties__(self.Fluid_hot)
        self.__update_fluid_properties__(self.Fluid_hot_liquid)

        # Cold fluid
        self.Fluid_cold         = Fluid(FluidsList.Air)
        self.Fluid_cold.update(Input.temperature(self.T_cold),Input.pressure(self.P_cold))
        self.__update_fluid_properties__(self.Fluid_cold)

    def __update_fluid_properties__(self,fluid):

        fluid.rho            = fluid.density
        fluid.Pr             = fluid.prandtl
        fluid.Re_D           = (fluid.density*self.u*self.outer_diameter)/self.nu
        fluid.k              = fluid.conductivity
        fluid.H              = fluid.enthalpy
        fluid.x              = fluid.quality

    def update_temperature(self,T):
        self.T = T
        self.fluid.update(Input.temperature(T))

    def __inner_convective_HTC__(self):
        h_hL           = (0.023*self.Fluid_hot_liquid.Re_D**(0.8)*self.Fluid_hot_liquid.Pr**(0.4)*self.Fluid_hot_liquid.k)/self.inner_diameter        
        self.h_h            = h_hL*((1-self.Fluid_hot.quality)**0.8+(3.8*self.Fluid_hot.quality**0.7 *(1-self.Fluid_hot.quality)**0.04)/(self.Fluid_hot.Pr**0.38))

    # For first row pipes
    def __outer_convective_HTC_row1__(self):
        Nu_D   = 0.3 + (0.62*self.Fluid_cold.Re_D**(1/2)*self.Fluid_cold.Pr**(1/3))/((1+(0.4/self.Fluid_cold.Pr)**(2/3))**(1/4)) * ((1+(self.Fluid_cold.Re_D)/282000)**(5/8))**(4/5)
        self.h_c    = (Nu_D*self.Fluid_cold.k)/self.outer_diameter

    # For all downstream pipes
    def __outer_convective_HTC_downstream__(self):



    def __inner_overall_HTC__(self):

        self.U_i = (1/(self.h_c) + 1/(self.h_h)  )**-1 

    def overall_HTC(self):
        self.__inner_convective_HTC__()
        self.__outer_convective_HTC__()
        self.__inner_overall_HTC__()

        self.U = 


Node.overall_HTC



    
    

    
    






# Define the fluid (water in this case)
water = Fluid(FluidsList.Water)

# Set the state using temperature (Celsius) and pressure (Pa)
water.update(Input.temperature(50), Input.pressure(1e5))  # 1 bar = 100000 Pa

# Get the specific enthalpy in J/kg
enthalpy = water.enthalpy

print(f"Enthalpy of water at 50°C and 1 bar: {enthalpy/1000:.2f} kJ/kg")