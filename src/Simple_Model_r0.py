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
    
    def with_state(self,Input1,Input2):
        '''Returns a new Fluid object with the specified state'''
        new_fluid = Fluid(self.name, self.m_dot)
        new_fluid.set_geometry(self.Geometry, self.flow_type)
        new_fluid.update(Input1, Input2)
        return new_fluid
    
    def liquid_phase(self):
        '''Returns a new Fluid object in the liquid phase at the current pressure'''
        if self.phase.name != "TwoPhase":
            raise ValueError("Fluid is not in a two-phase state.")
        
        liquid_phase = Fluid(self.name, self.m_dot)
        liquid_phase.set_geometry(self.Geometry, self.flow_type)
        liquid_phase.update(Input.quality(0), Input.pressure(self.pressure))
        return liquid_phase
    
    # Returns a new Fluid object in the vapor phase at the current pressure
    def vapor_phase(self):
        '''Returns a new Fluid object in the vapor phase at the current pressure'''
        if self.phase.name != "TwoPhase":
            raise ValueError("Fluid is not in a two-phase state.")
        
        vapor_phase = Fluid(self.name, self.m_dot)
        vapor_phase.update(Input.quality(1), Input.pressure(self.pressure))
        vapor_phase.set_geometry(self.Geometry,self.flow_type)
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

            A = (np.pi*self.Geometry.inner_diameter**2)/4
            u = self.m_dot / (self.density * A)
            Re = (self.density * self.Geometry.inner_diameter*u)/(self.dynamic_viscosity)
            
        elif self.flow_type == 'external':
            u = self.m_dot / (self.density * self.Geometry.duct_area)
        
            Re = u  * self.Geometry.outer_diameter / self.kinematic_viscosity
            
        else :
            raise ValueError("Invalid flow type. Use 'internal' or 'external'.")
        return Re
        
    # Add alias methods for commonly used properties
    # For creating alias properties
    def alias_property(original_name):
        return property(lambda self: getattr(self, original_name),  # Getter
                        lambda self, value: setattr(self, original_name, value) # Setter
                        )

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
    def __init__(self,Geometry,x_pos,node_length,T_hot_init,T_cold_init,P_hot_init,P_cold_init,m_dot_hot,m_dot_cold,is_boundary=False):
        self.Geometry = Geometry
        self.is_boundary = is_boundary
        self.x_pos = x_pos
        self.node_length = node_length

        # Hot fluid
        self.Fluid_hot = Fluid(FluidsList.Water, m_dot_hot)
        self.Fluid_hot.update(Input.temperature(T_hot_init), Input.pressure(P_hot_init))
        self.Fluid_hot.set_geometry(self.Geometry,type='internal')

        # Cold fluid
        self.Fluid_cold = Fluid(FluidsList.Air, m_dot_cold)
        self.Fluid_cold.update(Input.temperature(T_cold_init),Input.pressure(P_cold_init))
        self.Fluid_cold.set_geometry(self.Geometry, type='external')

    def __inner_convective_HTC__(self):
        
        # If not two-phase
        if self.Fluid_hot.phase.name != "TwoPhase":
            NotImplementedError("Single phase heat transfer not implemented yet.")
            self.h_h = -np.ln((-Tmx + T_s)/(-Tm_in + T_s))*m_dot*c_P/(P*x)
        else:
           # Shah correlation for boiling inside tubes
           liquid_phase = self.Fluid_hot.liquid_phase()

           # Shah 1979 eq 5
           h_L = (0.023*liquid_phase.Re**(0.8)*liquid_phase.Pr**(0.4)*liquid_phase.k)/self.Geometry.inner_diameter
           
           x = self.Fluid_hot.quality
           
           # Shah 1979 eq 8
           h_TP = h_L * ((1-x)**0.8 + (3.8*x**0.76 * (1-x)**0.04)/(self.Fluid_hot.Pr**0.38))
           self.h_h = h_TP
           
    def __outer_convective_HTC__(self):        
        Nu_D = 0.3 + (
            0.62*self.Fluid_cold.Re**(1/2)*self.Fluid_cold.Pr**(1/3))/((1+(0.4/self.Fluid_cold.Pr)**(2/3))**(1/4)) * ((1+(self.Fluid_cold.Re)/282000)**(5/8))**(4/5)
        self.h_c    = (Nu_D*self.Fluid_cold.k)/self.Geometry.outer_diameter

    def overall_HTC(self):
        self.__inner_convective_HTC__()
        self.__outer_convective_HTC__()
        
        # Inner overall heat transfer coefficient 
        # From heat transfer in single and multiphase systems (10.69) p. 495
        U_i = np.abs(np.pow((1/self.h_c + self.Geometry.inner_diameter/(self.Geometry.outer_diameter*self.h_h)),-1))
            
        inner_area = self.Geometry.inner_perimeter * self.node_length
        Delta_H =( U_i * inner_area )/self.Fluid_cold.m_dot * (self.Fluid_hot.temperature - self.Fluid_cold.temperature)
        
        return Delta_H
        
        
    def pressure_drop(self):
        pass
    
    
if __name__ == "__main__":
    geom = Geometry(
    # Height and length of the heat exchanger
    width=1,
    length=1,
    # Pipe dimensions
    pipe_outer_diameter=0.0603, 
    pipe_wall_thickness=0.005,
    # Arrengement dimensions
    arrangement="Inline", 
    transverse_pitch=0.04, 
    longitudinal_pitch=0.04
    )

    node1 = Node(
        Geometry=geom,
        x_pos=0,
        node_length=0.1, #m
        T_hot_init=343, # C 
        T_cold_init=200, # C 
        P_hot_init=158e5, # Pa
        P_cold_init=1e5, # Pa
        m_dot_hot=6, # kg/s
        m_dot_cold=50, # kg/s
        is_boundary=True
    )


    node1.Fluid_hot.update(Input.pressure(158e5), Input.enthalpy(node1.Fluid_hot.dew_point_at_pressure(158e5).enthalpy-1000))
    
    Delta_H = node1.overall_HTC()
    print(Delta_H)
    