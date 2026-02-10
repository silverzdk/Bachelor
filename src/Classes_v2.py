import numpy as np 
import scipy as sp
import matplotlib.pyplot as plt
from enum import StrEnum
from functools import cached_property 


import pyfluids as pf
from pyfluids import FluidsList, Input
# Set the units system to SI with Celsius
pf.PyFluidsConfig.units_system = pf.UnitsSystem.SIWithCelsius

g = 9.81 # m/s^2

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
            self.inner_diameter  = self.outer_diameter - 2 * self.wall_thickness
            self.inner_perimeter = np.pi * self.inner_diameter
            self.outer_perimeter = np.pi * self.outer_diameter
            
            self.inner_area = np.pi * (self.inner_diameter/2)**2


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

class Node: # Fluid state at a specific position. 
    def __init__(self,fluid_type,x_pos,mass_flow_rate,Input_1,Input_2):
        self.x_pos          = x_pos # m
        self.m_dot          = mass_flow_rate # kg/s

        self.fluid = pf.Fluid(fluid_type)
        self.fluid.update(Input_1, Input_2)

    def update_state(self, Input_1, Input_2):
        """Update thermodynamic state in-place."""
        self.fluid.update(Input_1, Input_2)

    def liquid_phase(self):
        """Return a Fluid object representing the saturated liquid phase."""
        if self.fluid.phase.name != "TwoPhase":
            raise ValueError("Node is not in a two-phase state")

        liquid = pf.Fluid(self.fluid.name)
        liquid.update(
            Input.pressure(self.fluid.pressure),
            Input.quality(0.0)
        )
        return liquid

    def vapor_phase(self):
        """Return a Fluid object representing the saturated vapor phase."""
        if self.fluid.phase.name != "TwoPhase":
            raise ValueError("Node is not in a two-phase state")

        vapor = pf.Fluid(self.fluid.name)
        vapor.update(
            Input.pressure(self.fluid.pressure),
            Input.quality(1.0)
        )
        return vapor

    def __repr__(self):
        return f"Fluid(type={self.fluid.name}, T={self.fluid.temperature} C, P={self.fluid.pressure} Pa, , H={self.fluid.enthalpy} J/kg, x={self.fluid.quality})"

class Internal: # Internal correlations
    def __init__(self,mass_flow,Node_in,Node_out,Geometry,Delta_x):
        self.Node_in    = Node_in
        self.Node_out   = Node_out
        self.Geometry   = Geometry
        self.m_dot      = mass_flow
        self.L          = Delta_x

    def pressure_drop(self):
        '''Determines the type of pressure drop based on the flow type and phase change'''
        phase_in    = self.Node_in.fluid.phase.name
        phase_out   = self.Node_out.fluid.phase.name

        if phase_in in ("Liquid", "Gas") and phase_out in ("Liquid", "Gas"):
            dp = self.monophase_pressure_drop(self.L,self.Node_in,self.Node_out)

        elif phase_in == "TwoPhase" and phase_out == "TwoPhase":
            dp = self.two_phase_pressure_drop(self.Node_in.fluid.quality,self.Node_out.fluid.quality,self.Node_in.liquid_phase(),self.Node_out.liquid_phase(),self.Node_in.vapor_phase(),self.Node_out.vapor_phase())

        elif phase_in != phase_out:
            
            # Phase change direction
            # Representative pressure drop for phase change segment
            p_avg = (self.Node_in.fluid.pressure + self.Node_out.fluid.pressure) / 2

            if phase_in == "Gas" and phase_out == "TwoPhase":
                
                sat_vap = Node(
                    fluid_type=self.Node_in.fluid.name,
                    x_pos=self.Node_in.x_pos,
                    mass_flow_rate=self.Node_in.m_dot,
                    Input_1=Input.quality(1.0),
                    Input_2=Input.pressure(p_avg)
                )
                h_sat_vap = sat_vap.fluid.enthalpy

                # Fraction of CV that is gas
                f = (self.Node_in.fluid.enthalpy - h_sat_vap) / (
                        self.Node_in.fluid.enthalpy - self.Node_out.fluid.enthalpy
                    )
                f = min(max(f, 0.0), 1.0)

                L_Gas = self.L * f
                dp_Gas = self.monophase_pressure_drop(
                    L_Gas,
                    self.Node_in,
                    sat_vap
                )

                dp_TwoPhase = self.two_phase_pressure_drop(
                    sat_vap.fluid.quality,
                    self.Node_out.fluid.quality,
                    sat_vap.liquid_phase(),
                    self.Node_out.liquid_phase(),
                    sat_vap.vapor_phase(),
                    self.Node_out.vapor_phase()
                )

                dp = dp_Gas + dp_TwoPhase

            elif phase_in == "TwoPhase" and phase_out == "Liquid":

                # Saturation at representative pressure (saturated liquid)
                sat_liq = Node(
                    fluid_type=self.Node_in.fluid.name,
                    x_pos=self.Node_in.x_pos,
                    mass_flow_rate=self.Node_in.m_dot,
                    Input_1=Input.quality(0.0),
                    Input_2=Input.pressure(p_avg)
                )
                h_sat_liq = sat_liq.fluid.enthalpy

                # Fraction of CV that is two-phase
                f = (self.Node_in.fluid.enthalpy - h_sat_liq) / (
                        self.Node_in.fluid.enthalpy - self.Node_out.fluid.enthalpy
                    )
                f = min(max(f, 0.0), 1.0)

                # Lengths of sub-sections
                L_Liquid   = self.L * (1.0 - f)

                # Two-phase pressure drop (from inlet quality → 0)
                dp_TwoPhase = self.two_phase_pressure_drop(
                    self.Node_in.fluid.quality,
                    sat_liq.fluid.quality,   # = 0.0
                    self.Node_in.liquid_phase(),
                    sat_liq.liquid_phase(),
                    self.Node_in.vapor_phase(),
                    sat_liq.vapor_phase()
                )

                # Liquid pressure drop (from saturation → outlet)
                dp_Liquid = self.monophase_pressure_drop(
                    L_Liquid,
                    sat_liq,
                    self.Node_out
                )

                dp = dp_TwoPhase + dp_Liquid

        return dp
    
    def monophase_pressure_drop(self,Length,Node_in,Node_out):
        print("Monophase check")

        # Calculate flow velocities
        u_in       = self.m_dot / (Node_in.fluid.density * self.Geometry.Tube.inner_area)
        u_out      = self.m_dot / (Node_out.fluid.density * self.Geometry.Tube.inner_area)

        # Calculate average properties
        rho_avg    = (Node_in.fluid.density + Node_out.fluid.density) / 2
        u_avg      = self.m_dot / (rho_avg * self.Geometry.Tube.inner_area)
        mu_avg     = (Node_in.fluid.dynamic_viscosity + Node_out.fluid.dynamic_viscosity) / 2           
        Re_avg     = self.Reynolds_number(rho_avg,u_avg,mu_avg,self.Geometry.Tube.inner_diameter)

        # Acceleration pressure drop
        DP_A        = 0.5 * (Node_out.fluid.density * u_out**2 - Node_in.fluid.density * u_in**2)
        
        # Frictional pressure drop (Darcy-Weisbach equation)
        # Calculate friction factor depending on flow regime
        zeta_fr = self.friction_factor(Re_avg)
        
        # Calculate frictional pressure drop
        DP_R    = zeta_fr * (Length/self.Geometry.Tube.inner_diameter) * 0.5 * (rho_avg * u_avg**2)

        # Calculate gravitational pressure drop
        DP_G    = 0 # Assuming horizontal flow for now

        # Total pressure drop of internal segment
        DP_tot  = DP_A + DP_R + DP_G

        return(DP_tot)    
    
    def two_phase_pressure_drop(self,x_in,x_out,Node_in_liq,Node_out_liq,Node_in_vap,Node_out_vap):
        print("two phase check 1")

        # Calculate average flow velocities
        u_liq_avg   = self.m_dot / (( (Node_in_liq.density + Node_out_liq.density) / 2) * self.Geometry.Tube.inner_area)
        u_vap_avg   = self.m_dot / (( (Node_in_vap.density + Node_out_vap.density) / 2) * self.Geometry.Tube.inner_area)

        # Calculate average properties
        rho_avg_liq = (Node_in_liq.density + Node_out_liq.density) / 2
        rho_avg_vap = (Node_in_vap.density + Node_out_vap.density) / 2
        mu_avg_liq  = (Node_in_liq.dynamic_viscosity + Node_out_liq.dynamic_viscosity) / 2
        mu_avg_vap  = (Node_in_vap.dynamic_viscosity + Node_out_vap.dynamic_viscosity) / 2

        # Calculate Reynolds numbers for liquid and vapor phases
        Re_avg_liq  = self.Reynolds_number(rho_avg_liq,u_liq_avg,mu_avg_liq,self.Geometry.Tube.inner_diameter)
        Re_avg_vap  = self.Reynolds_number(rho_avg_vap,u_vap_avg,mu_avg_vap,self.Geometry.Tube.inner_diameter)

        # Friction factor for liquid and vapor phases
        zeta_fr_liq = self.friction_factor(Re_avg_liq)
        zeta_fr_vap = self.friction_factor(Re_avg_vap)

        M_dot = self.m_dot / self.Geometry.Tube.inner_area

        # Pressure drop of liquid and vapor phases
        DP_R_liq    = zeta_fr_liq * (M_dot**2)/(2 * rho_avg_liq * self.Geometry.Tube.inner_diameter)
        DP_R_vap    = zeta_fr_vap * (M_dot**2)/(2 * rho_avg_vap * self.Geometry.Tube.inner_diameter)

        # Assigning coefficients
        A = DP_R_liq
        B = DP_R_vap

        def F(A,B,x):
            F_ABx = -(3/4) * (1-x)**(3/4) * (A+2*(B-A)*x)+(1/4)*B*x**4 -9/14 * (B-A)*(1-x)**(7/3)
            return F_ABx
        
        DP_R = F(A,B,x_out) - F(A,B,x_in)

        return(DP_R)

    def Reynolds_number(self,rho,u,mu,d):
        Re = (rho * u * d) / mu
        return Re

    def friction_factor(self,Re):
        if 3000 < Re < 100000:
            zeta = 0.3614/((Re**(1/4)))
        elif 2*10**4 < Re < 2*10**6:
            zeta = 0.3614/((Re**(1/4)))
        elif Re > 10**6:
            f = lambda z: 1/z + 0.8 - 2*np.log(Re*np.sqrt(z))
            zeta = sp.optimize.fsolve(f, x0=0.001)[0]
        return zeta

    def heat_transfer(self):
        pass

class Segment: # External correlations
    def __init__(self,rownumber,node,geometry):
        self.rownumber = rownumber # Row number in the bank
        self.node = node # Node object representing the fluid state
        self.geometry = geometry # Geometry class

    def Flow_velocity(self):
        '''Calculates the flow velocity in the cold segment'''
        u_g = self.node.m_dot / (self.geometry.Free_flow_area * self.node.pressure)
        return u_g




