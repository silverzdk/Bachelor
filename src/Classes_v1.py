import numpy as np 
import scipy as sp
import matplotlib.pyplot as plt
from enum import StrEnum
from functools import cached_property 


import pyfluids as pf
from pyfluids import FluidsList, Input
# Set the units system to SI with Celsius
pf.PyFluidsConfig.units_system = pf.UnitsSystem.SIWithCelsius

# Enum classes for various types
class FlowType(StrEnum):
    '''Enum for flow types'''
    Internal = "internal"
    External = "external"
    
class arrangementType(StrEnum):
    '''Enum for tube arrangement types'''
    Inline      = "Inline"
    Staggered   = "Staggered"
    Single_Tube = "Single_Tube"
    Single_Row  = "Single_Row"
    
class finType(StrEnum):
    '''Enum for fin types'''
    Rectangular = "rectangular"
    Circular    = "circular"
    Spiral      = "spiral"

# Classes for heat exchanger modeling

# Fluid class that extends pyfluids.Fluid
class Fluid(pf.Fluid):
    def __init__(self,fluid_type,mass_flow_rate=None):
        # Initialize the base Fluid class
        super().__init__(fluid_type)
        self.Geometry = None
        self.mass_flow_rate = mass_flow_rate

    def set_geometry(self,Geometry,type=FlowType.Internal):
        '''Sets the geometry of the heat exchanger for flow calculations
        args:
            Geometry: Geometry object containing geometry information
            type: 'internal' or 'external' flow type
        '''
        self.Geometry = Geometry
        self.flow_type = type
        
        if type == FlowType.Internal:
            self.characteristic_length = Geometry.inner_diameter
        if type == FlowType.External:
            self.characteristic_length = Geometry.outer_diameter

    def with_state(self,Input1,Input2):
        '''Returns a new Fluid object with the specified state'''
        if self.m_dot is not None:
            new_fluid = Fluid(self.name, self.m_dot)
        else: 
            new_fluid = Fluid(self.name)
        if self.Geometry is not None:
            new_fluid.set_geometry(self.Geometry, self.flow_type)
        
        new_fluid.update(Input1, Input2)
        return new_fluid

    def liquid_phase(self):
        '''Returns a new Fluid object in the liquid phase at the current pressure'''
        if self.phase.name != "TwoPhase":
            raise ValueError(f"{self} is not in a two-phase state.") 

        liquid_phase = self.with_state(Input.quality(0), Input.pressure(self.pressure))
       
        return liquid_phase
    
    # Returns a new Fluid object in the vapor phase at the current pressure
    def vapor_phase(self):
        '''Returns a new Fluid object in the vapor phase at the current pressure'''
        
        if self.phase.name != "TwoPhase":
            raise ValueError(f"{self} is not in a two-phase state.") 
        
        vapor_phase = self.with_state(Input.quality(1), Input.pressure(self.pressure))

        return vapor_phase
    
    def __repr__(self):
        return f"Fluid(type={self.name}, T={self.temperature} C, P={self.pressure} Pa, , H={self.enthalpy} J/kg, x={self.quality})"
    
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
    class _Tube:
        '''Inner class representing tube geometry'''
        def __init__(self):
            self.geom_set = False
        def set_geometry(self,outer_diameter,wall_thickness):
            '''Sets the geometry of the tube
            args:   
                outer_diameter: Outer diameter of the tube (m)
                wall_thickness: Wall thickness of the tube (m)
            '''
            self.geom_set = True
            self.outer_diameter = outer_diameter
            self.wall_thickness = wall_thickness
            
            # Derived properties
            self.inner_diameter = self.outer_diameter - 2 * self.wall_thickness
            self.inner_perimeter= np.pi * self.inner_diameter
            self.outer_perimeter= np.pi * self.outer_diameter
            
    class _Duct:
        '''Inner class representing duct geometry'''
        def __init__(self):
            self.geom_set = False
        def set_geometry(self,a,b):
            '''Sets the geometry of the duct
            args:
                a: Duct width (m)
                b: Duct height (m)
            '''
            self.geom_set = True
            self.a = a
            self.b  = b
            
            # Derived properties
            self.area   = self.a * self.b
    
    class _Fin:
        '''Inner class representing fin geometry'''
        def __init__(self):
            self.geom_set = False
        def set_geometry(self,average_fin_thickness,fin_height,fin_spacing,fin_type=finType.Rectangular):
            '''Sets the geometry of the fin
            args:
                average_fin_thickness: Average thickness of the fin (m)
                fin_height: Height of the fin measured from tube OD (m)
                fin_spacing: Spacing between fins / fin pitch (m)
                fin_type: Type of the fin from finType enum (default is finType.Rectangular)
            '''
            self.geom_set = True
            self.average_fin_thickness  = average_fin_thickness #delta_r
            self.fin_height             = fin_height #l_r
            self.fin_spacing            = fin_spacing #s_r
            self.fin_type               = fin_type
            
            match fin_type:
                case finType.Rectangular:
                    self.square_height = self.fin_height*2 + self.Tube.outer_diameter
                case finType.Circular:
                    self.finning_diameter = self.Fin.fin_height*2 + self.Tube.outer_diameter
                case finType.Spiral:
                    raise NotImplementedError("Spiral fin type not implemented yet.")
    
    class _Bank: 
        '''Inner class representing bank geometry'''
        def __init__(self):
            self.geom_set = False
        def set_geometry(self, transverse_number_of_rows,longitudinal_number_of_rows,transverse_pitch,longitudinal_pitch,arrangement=arrangementType.Inline):
            ''''Sets the geometry of the tube bank
            args:
                transverse_number_of_rows: Number of tube rows in the bank transverse to flow direction (-)
                longitudinal_number_of_rows: Number of tube rows in the bank longitudinal to flow direction (-)
                transverse_pitch: Transverse pitch of tubes (m)
                longitudinal_pitch: Longitudinal pitch of tubes (m)
                arrangement: Arrangement type of tubes (Inline, Staggered, Single_Tube, Single_Row) (-)
            '''
            self.geom_set = True
            self.transverse_number_of_rows     = transverse_number_of_rows
            self.longitudinal_number_of_rows   = longitudinal_number_of_rows
            self.transverse_pitch   = transverse_pitch
            self.longitudinal_pitch = longitudinal_pitch
            self.arrangement        = arrangement
            
            self.total_number_of_tubes = self.transverse_number_of_rows * self.longitudinal_number_of_rows 
            
            # Derived properties
            match arrangement:
                case arrangementType.Inline:
                    self.diagonal_pitch = None
                case arrangementType.Staggered:
                    self.diagonal_pitch = np.sqrt(self.longitudinal_pitch**2+(self.transverse_pitch/2)**2)
                case arrangementType.Single_Tube:
                    self.diagonal_pitch = None
                case arrangementType.Single_Row:
                    self.diagonal_pitch = None
    
    def __init__(self):
        # Initialize inner classes
        self.Tube = self._Tube()
        self.Duct = self._Duct()
        self.Fin  = self._Fin()
        self.Bank = self._Bank()
    
    def minimum_free_flow_area(self):
        '''Calculates the minimum free flow area based on the geometry parameters
            Returns:
            F: Minimum free flow area (m^2)'''
        
        # Check that all necessary geometry parameters are set
        if not (self.Tube.geom_set and self.Duct.geom_set and self.Fin.geom_set and self.Bank.geom_set):
            raise ValueError("All geometry parameters must be set before calculating minimum free flow area.")
        
        # Calculate equivalent diameter
        d_cl = self.Tube.outer_diameter + (2* self.Fin.fin_height*self.Fin.average_fin_thickness) / self.Fin.fin_spacing
        
        S_1 = self.Bank.transverse_pitch
        S_2 = self.Bank.longitudinal_pitch

        # ϕ parameter:
        phi_cl = (S_1 - d_cl) / (S_2 - d_cl)

        # Free flow area:
        # Note currently uses duct b as tube length - may need to be changed
        if self.Bank.arrangement == arrangementType.Inline or phi_cl <= 2:
            F = self.Duct.area - self.Bank.transverse_number_of_rows * d_cl * self.Duct.b
        elif phi_cl > 2:
            F = (self.Duct.area - self.Bank.transverse_number_of_rows * d_cl * self.Duct.b) * 2/ phi_cl

        return F
    
    @cached_property
    def sigma_1(self):
        '''Calculates the transverse pitch to diameter ratio'''
        return self.Bank.transverse_pitch / self.Tube.outer_diameter
    
    @cached_property
    def sigma_2(self):
        '''Calculates the longitudinal pitch to diameter ratio'''
        return self.Bank.longitudinal_pitch / self.Tube.outer_diameter
    
    @cached_property
    def Psi_r(self):
        '''Calculates the fin coefficient for finned tubes'''
        # In the below calculations we assume constant fin thickness so that the average fin thickness equals the fin thickness at base and tip
        match self.Fin.fin_type:
            case finType.Rectangular:
                Psi_r = (2*(self.Fin.square_height**2 - 0.785 * self.Tube.outer_diameter**2 + 2 * self.Fin.square_height*self.Fin.average_fin_thickness))(np.pi * self.Tube.outer_diameter * self.Fin.fin_spacing) + (1 - self.Fin.average_fin_thickness/self.Fin.fin_spacing)
            case finType.Circular: 
                Psi_r = 1/(2*self.Tube.outer_diameter*self.Fin.fin_spacing) *(self.Fin.finning_diameter**2 - self.Tube.outer_diameter**2 + 2*self.Fin.finning_diameter * self.Fin.average_fin_thickness) + (1- self.Fin.average_fin_thickness/self.Fin.fin_spacing)
    
    @cached_property
    def A_r(self):
        '''Calculates the heat transfer area for fins'''
        
        match self.Fin.fin_type:
            case finType.Rectangular:
                A_r = 2 * (self.Fin.square_height**2 - 0.785 * self.Tube.outer_diameter**2 + 2 * self.Fin.square_height*self.Fin.average_fin_thickness) * self.Bank.finned_tube_segment / self.Fin.fin_spacing * self.Bank.total_number_of_tubes
                return A_r
            case finType.Spiral | finType.Circular:
                A_r = np.pi/2 * (self.Fin.finning_diameter**2 - self.Tube.outer_diameter**2 + 2 * self.Fin.finning_diameter * self.Fin.average_fin_thickness) * self.Bank.finned_tube_segment / self.Fin.fin_spacing * self.Bank.total_number_of_tubes
        
        
    
    @cached_property
    def A_t(self):
        '''Calculates the heat transfer area for tubes'''
        raise NotImplementedError("Tube area calculation not implemented yet.")
    
    @cached_property
    def A(self):
        '''Calculates the total heat transfer area'''
        return self.A_r + self.A_t
    
class Node:

    def __init__(self,Geometry,x_pos,node_length,T_hot_init,T_cold_init,P_hot_init,P_cold_init,m_dot_hot,m_dot_cold,is_boundary=False):
        self.Geometry       = Geometry
        self.is_boundary    = is_boundary
        self.x_pos          = x_pos
        self.node_length    = node_length

        # Hot fluid
        self.Fluid_hot = Fluid(FluidsList.Water, m_dot_hot)
        self.Fluid_hot.update(Input.temperature(T_hot_init), Input.pressure(P_hot_init))
        self.Fluid_hot.set_geometry(self.Geometry,type='internal')

        # Cold fluid
        self.Fluid_cold = Fluid(FluidsList.Air, m_dot_cold)
        self.Fluid_cold.update(Input.temperature(T_cold_init),Input.pressure(P_cold_init))
        self.Fluid_cold.set_geometry(self.Geometry, type='external')
        
        # Unitialized variables 
        self.h_c = None
        self.h_h = None

    def __Re_eq__():
        
        #Conventional diameter of finned tube:
        d_cl = d + (2*l_r*delta_r)/(s_r)

        #phi parameter:
        phi_cl = (S_1-d_cl)/(S_2-d_cl)

        #Free flow area:
        if phi_cl <= 2:
            F = a*b-z_1*L_ccrs*d_cl
        elif phi_cl > 2:
            F = (a*b - z_1*L_ccrs*d_cl)*(2/phi_cl)

        #Design gas velocity:
        u_g     = G_g/(rho_g * F)

        #Reynold's number is calculated:
        Re_eq = (u_g*d_eq)/v_g

    def __inner_convective_HTC__(self):
        
        # If not two-phase
        if self.Fluid_hot.phase.name != "TwoPhase":
            NotImplementedError("Single phase heat transfer not implemented yet.")
            
        #     self.h_h = -np.ln((-Tmx + T_s)/(-T + T_s))*m_dot*c_P/(P*x)

        # Shah correlation for boiling inside tubes
        liquid_phase = self.Fluid_hot.liquid_phase()

        # Shah 1979 eq 5
        h_L = (0.023*liquid_phase.Re**(0.8)*liquid_phase.Pr**(0.4)*liquid_phase.k)/self.Geometry.inner_diameter
           
        x = self.Fluid_hot.quality
        p_crit   = self.Fluid_hot.critical_pressure
        p_r      = self.Fluid_hot.pressure / p_crit

        # Shah 1979 eq 8
        h_TP = h_L * ((1-x)**0.8 + (3.8*x**0.76 * (1-x)**0.04)/(p_r**0.38))
        self.h_h = h_TP
        
    def __outer_convective_HTC__(self):
        
        F = self.Geometry.minimum_free_flow_area()
        u_g = self.Fluid_cold.m_dot / (self.Fluid_cold.rho * F)
        
        match self.Geometry.Bank.arrangement:
            case arrangementType.Inline:
                X = 4 * (2 + self.Geometry.Psi_r/7 - self.Geometry.sigma_2) 
            case arrangementType.Staggered:
                X = self.Geometry.sigma_1 / self.Geometry.sigma_2  - 1.26 / self.Geometry.Psi_r - 2
            case _: 
                NotImplementedError("Only Inline and Staggered tube arrangements are implemented for external heat transfer coefficient calculation.")
        
        n = 0.7 + 0.08 * np.tanh(X) + 0.005 * self.Geometry.Psi_r
        C_q = (1.36-np.tanh(X))*( (1.1)/(self.Geometry.Psi_r+8) - 0.014)
        
        # Calculate C_z based on number of rows and arrangement
        if self.Geometry.Bank.transverse_number_of_rows >= 8:
            C_z = 1
        elif self.Geometry.Bank.arrangement == arrangementType.Inline or (self.Geometry.sigma_1/self.Geometry.sigma_2) >= 2:
            C_z = 3.5 * self.Geometry.Bank.transverse_number_of_rows**(0.03) - 2.72
        elif self.Geometry.Bank.arrangement == arrangementType.Staggered and (self.Geometry.sigma_1/self.Geometry.sigma_2) < 2:
            C_z = 3.15 * self.Geometry.Bank.transverse_number_of_rows**(0.05) - 2.5
        else:
            raise ValueError("Error in calculating C_z for external heat transfer coefficient.")
        
        h_c = 1.13 * C_z * C_q * (self.Fluid_cold.k / self.Geometry.outer_diameter) * (u_g * self.Geometry.outer_diameter / self.Fluid_cold.nu)**n * self.Fluid_cold.Pr**0.33
        
        self.h_c = h_c
        
    # def __Reduced_HTC__(self):

    #     if self.h_c is None:
    #         self.__outer_convective_HTC__()
            
    #     h_1rdc = (self.Geometry)
        
        
        


class Row:

    def __init__(self,Geometry,longitudinal_row_pos,node_length,T_hot_init,T_cold_init,P_hot_init,P_cold_init,m_dot_hot,m_dot_cold,roughness,is_boundary=False):
        self.Geometry       = Geometry
        self.is_boundary    = is_boundary
        self.longitudinal_row_pos        = longitudinal_row_pos


        # Hot fluid
        self.Fluid_hot = Fluid(FluidsList.Water, m_dot_hot)
        self.Fluid_hot.update(Input.temperature(T_hot_init), Input.pressure(P_hot_init))
        self.Fluid_hot.set_geometry(self.Geometry,type='internal')

        # Cold fluid
        self.Fluid_cold = Fluid(FluidsList.Air, m_dot_cold)
        self.Fluid_cold.update(Input.temperature(T_cold_init),Input.pressure(P_cold_init))
        self.Fluid_cold.set_geometry(self.Geometry, type='external')


if __name__ == "__main__":
    Geom = Geometry()
    Geom.minimum_free_flow_area()