import numpy as np 
import scipy as sp
import matplotlib.pyplot as plt

import pyfluids as pf
from pyfluids import FluidsList, Input


# Fluid class that extends pyfluids.Fluid
class Fluid(pf.Fluid):
    def __init__(self,fluid_type,mass_flow_rate=None):
        # Initialize the base Fluid class
        super().__init__(fluid_type)
        self.Geometry = None
        self.mass_flow_rate = mass_flow_rate
    
    def set_geometry(self,Geometry,type='internal'):
        '''Sets the geometry of the heat exchanger for flow calculations
        args:
            Geometry: Geometry object containing geometry information
            type: 'internal' or 'external' flow type
        '''
        self.Geometry = Geometry
        self.flow_type = type
        
        if type == 'internal':
            self.characteristic_length = Geometry.inner_diameter
        if type == 'external':
            self.characteristic_length = Geometry.outer_diameter

    def liquid_phase(self):
        '''Returns a new Fluid object in the liquid phase at the current pressure'''
        if self.phase.name != "TwoPhase":
            raise ValueError("Fluid is not in a two-phase state.")
        
        liquid_phase = Fluid(self.name, self.m_dot)
        liquid_phase.update(pf.Input.quality(0), pf.Input.pressure(self.pressure))
        liquid_phase.set_geometry(self.Geometry)
        return liquid_phase
    
    # Returns a new Fluid object in the vapor phase at the current pressure
    def vapor_phase(self):
        '''Returns a new Fluid object in the vapor phase at the current pressure'''
        if self.phase.name != "TwoPhase":
            raise ValueError("Fluid is not in a two-phase state.")
        
        vapor_phase = Fluid(self.name, self.m_dot)
        vapor_phase.update(pf.Input.quality(1), pf.Input.pressure(self.pressure))
        vapor_phase.set_geometry(self.Geometry)
        return vapor_phase
    
    def __repr__(self):
        return f"Fluid(type={self.name}, T={self.temperature} C, P={self.pressure} Pa, , H={self.enthalpy} J/kg, x={self.quality})"
    
    # Define additional properties for convenience
    @property
    def Re(self):
        '''Returns the Reynolds number based on the current mass flow rate and pipe geometry'''
        if self.Geometry is None:
            raise ValueError("Pipe geometry not set. Use set_geometry() method.")
        if self.m_dot is None:
            raise ValueError("Mass flow rate not set.")
        
        # Calculate based on flow type
        if self.flow_type == 'internal':
            
            u = self.m_dot / (self.density * (np.pi * (self.Geometry.inner_diameter / 2)**2))
            Re = (self.density * u * self.Geometry.inner_diameter) / self.dynamic_viscosity
            
        elif self.flow_type == 'external':
            V = self.m_dot / (self.density * self.Geometry.duct_area)
            
            if self.Geometry.arrangement == "Inline":
                V_max = V* (self.Geometry.transverse_pitch / (self.Geometry.transverse_pitch - self.Geometry.outer_diameter))
            elif self.Geometry.arrangement == "Staggered":
                V_max = V* (self.Geometry.transverse_pitch/ (self.Geometry.diagonal_pitch - self.Geometry.outer_diameter))
        
        #-- UNFINISHED EXTERNAL FLOW CALCULATION --#
            
        else :
            raise ValueError("Invalid flow type. Use 'internal' or 'external'.")
        return Re
    
    @property
    def nu(self):
        return self.dynamic_viscosity / self.density
    
    # Add alias methods for commonly used properties
    # For creating alias properties
    def alias_property(original_name):
        return property(lambda self: getattr(self, original_name))

    k = alias_property("conductivity")
    H = alias_property("enthalpy")
    x = alias_property("quality")
    rho = alias_property("density")
    Pr = alias_property("prandtl")
    mu = alias_property("dynamic_viscosity")
    m_dot = alias_property("mass_flow_rate")

class Geometry:
    '''Class representing the geometry of the heat exchanger with given dimensions'''
    def __init__(self,width,length,pipe_outer_diameter,pipe_wall_thickness,transverse_pitch,longitudinal_pitch,arrangement="Inline",angle=0):
        self.width = width
        self.length = length
        self.outer_diameter = pipe_outer_diameter
        self.wall_thickness = pipe_wall_thickness
        self.transverse_pitch = transverse_pitch
        self.longitudinal_pitch = longitudinal_pitch
        self.arrangement = arrangement
        self.angle          = angle
        
        # Calculate derived properties
        # tube properties
        self.inner_diameter  = self.outer_diameter - 2 * self.wall_thickness
        self.inner_perimeter = np.pi * self.inner_diameter
        self.outer_perimeter = np.pi * self.outer_diameter
        
        # tube bank properties
        if arrangement == "Inline":
            self.diagonal_pitch = None
        elif arrangement == "Staggered":
            self.diagonal_pitch = np.sqrt(self.longitudinal_pitch**2+(self.transverse_pitch/2)**2)
        
        # Duct properties
        self.duct_area       = self.width * self.length

class Node:
    def __init__(self,pipe,x_pos,T_hot_init,T_cold_init,P_hot_init,P_cold_init,m_dot_hot,m_dot_cold,is_boundary=False):
        self.Pipe           = pipe
        self.is_boundary    = is_boundary
        self.x_pos          = x_pos

        # Hot fluid
        self.Fluid_hot = Fluid(FluidsList.Water, m_dot_hot)
        self.Fluid_hot.update(Input.temperature(T_hot_init), Input.pressure(P_hot_init))
        self.Fluid_hot.set_geometry(self.Pipe)

        # Cold fluid
        self.Fluid_cold = Fluid(FluidsList.Air, m_dot_cold)
        self.Fluid_cold.update(Input.temperature(T_cold_init),Input.pressure(P_cold_init))

    def __inner_convective_HTC__(self):
        
        # Shah correlation for boiling inside tubes
        liquid_phase = self.Fluid_hot.liquid_phase()
        
        h_hL = (0.023*liquid_phase.Re_D**(0.8)*liquid_phase.Pr**(0.4)*liquid_phase.k)/self.inner_diameter        
        self.h_h = h_hL*((1-self.Fluid_hot.quality)**0.8+(3.8*self.Fluid_hot.quality**0.7 *(1-self.Fluid_hot.quality)**0.04)/(self.Fluid_hot.Pr**0.38))

    # For first row pipes
    def __outer_convective_HTC_row1__(self):
        # Kapitel 7-4 i bogen s. 493 pdf Heat and mass transfer - Fundamentals and applications 
        Nu_D = 0.3 + (
            0.62*self.Fluid_cold.Re_D**(1/2)*self.Fluid_cold.Pr**(1/3))/((1+(0.4/self.Fluid_cold.Pr)**(2/3))**(1/4)) * ((1+(self.Fluid_cold.Re_D)/282000)**(5/8))**(4/5)
        self.h_c    = (Nu_D*self.Fluid_cold.k)/self.outer_diameter

    # For all downstream pipes
    def __outer_convective_HTC_downstream__(self):
        pass


    def __inner_overall_HTC__(self):

        self.U_i = (1/(self.h_c) + 1/(self.h_h)  )**-1 

    def overall_HTC(self):
        self.__inner_convective_HTC__()
        self.__outer_convective_HTC__()
        self.__inner_overall_HTC__()
