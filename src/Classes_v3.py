import numpy as np 
import scipy as sp
import matplotlib.pyplot as plt
from enum import StrEnum
from functools import cached_property 

import pyfluids as pf
from pyfluids import FluidsList, Input
# Set the units system to SI with Celsius
pf.PyFluidsConfig.units_system = pf.UnitsSystem.SIWithCelsius

# Global variables
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

# Geometry class with inner classes for tube, duct, fin, and bank geometry
class Geometry:
    """Class representing the geometry of a heat exchanger with defined dimensions."""

    class _Tube:
        """Inner class representing tube geometry."""

        def __init__(self,parent):
            self.parent = parent
            self.geom_set = False

        def set_geometry(self, outer_diameter, wall_thickness):
            """Set the geometric properties of the tube.

            Args:
                outer_diameter (float): Outer diameter of the tube in meters.
                wall_thickness (float): Wall thickness of the tube in meters.

            Notes:
                Derived properties such as inner diameter, perimeters, and
                internal area are computed automatically.
            """
            self.geom_set       = True
            self.outer_diameter = outer_diameter  
            self.wall_thickness = wall_thickness 

            # Derived properties
            self.inner_diameter     = self.outer_diameter - 2 * self.wall_thickness
            self.inner_perimeter    = np.pi * self.inner_diameter
            self.outer_perimeter    = np.pi * self.outer_diameter
            self.inner_area         = np.pi * (self.inner_diameter/2)**2

    class _Duct:
        """Inner class representing duct geometry."""

        def __init__(self,parent):
            self.parent = parent
            self.geom_set = False

        def set_geometry(self, a, b):
            """Set the geometric properties of the duct.

            Args:
                a (float): Duct height in meters.
                b (float): Duct width in meters.

            Notes:
                The frontal cross-sectional area is computed automatically.
            """
            self.geom_set = True
            self.a  = a
            self.b  = b

            # Derived properties
            self.area = self.a * self.b

    class _Fin:
        """Inner class representing fin geometry."""

        def __init__(self,parent):
            self.parent = parent
            self.geom_set = False

        def set_geometry(self, average_fin_thickness, fin_height, fin_spacing, fin_type=finType.Rectangular):
            """Set the geometric properties of the fins.

            Args:
                average_fin_thickness (float): Average fin thickness in meters.
                fin_height (float): Height of the fin measured from tube OD in meters.
                fin_spacing (float): Spacing between fins (fin pitch) in meters.
                fin_type (finType): Type of fin geometry. Defaults to finType.Rectangular.

            Raises:
                NotImplementedError: If the fin type is Spiral.

            Notes:
                Additional derived geometric parameters are computed depending
                on the selected fin type.
            """
            self.geom_set = True
            self.average_fin_thickness  = average_fin_thickness
            self.fin_height             = fin_height
            self.fin_spacing            = fin_spacing
            self.fin_type               = fin_type
        
        @cached_property
        def square_height(self):
            """Compute the square fin height.

            Returns:
                float: Tube outer diameter plus twice the fin height.
            Raises:
                ValueError: If the fin type is not Rectangular.
                ValueError: If the tube geometry is not set before computing square height.

            """
            if not self.fin_type == finType.Rectangular:
                raise ValueError("Square height is only defined for rectangular fins.")
            if not self.parent.Tube.geom_set:
                raise ValueError("Tube geometry must be set before computing square height.")
            
            return self.parent.Tube.outer_diameter+2*self.fin_height
        
        
        @cached_property
        def finning_diameter(self):
            """Compute the finning diameter.

            Returns:
                float: Tube outer diameter plus twice the fin height.
            Raises:
                ValueError: If the fin type is not Circular.
                ValueError: If the tube geometry is not set before computing finning diameter.

            """
            if not self.fin_type == finType.Circular:
                raise ValueError("Finning diameter is only defined for circular fins.")
            if not self.parent.Tube.geom_set:
                raise ValueError("Tube geometry must be set before computing finning diameter.")
            
            return self.parent.Tube.outer_diameter+2*self.fin_height
        
        @cached_property
        def Psi_r(self):
            """Compute the fin efficiency coefficient for finned tubes.

            Returns:
                float: Fin coefficient Psi_r.
                
            Raises:
                ValueError: If the tube geometry is not set before computing Psi_r.

            Notes:
                Assumes constant fin thickness for simplification.
            """
            
            if not self.parent.Tube.geom_set:
                raise ValueError("Tube geometry must be set before computing Psi_r.")
            
            match self.fin_type:
                case finType.Rectangular:
                    Psi_r = (2*(self.square_height**2 - 0.785 * self.parent.Tube.outer_diameter**2
                            + 2 * self.square_height*self.average_fin_thickness)) \
                            / (np.pi * self.parent.Tube.outer_diameter * self.fin_spacing) \
                            + (1 - self.average_fin_thickness/self.fin_spacing)
                case finType.Circular: 
                    Psi_r = 1/(2*self.parent.Tube.outer_diameter*self.fin_spacing) * \
                            (self.finning_diameter**2 - self.parent.Tube.outer_diameter**2
                            + 2*self.finning_diameter * self.average_fin_thickness) \
                            + (1 - self.average_fin_thickness/self.fin_spacing)
            return Psi_r

    class _Bank:
        """Inner class representing bank geometry."""

        def __init__(self,parent):
            self.parent = parent
            self.geom_set = False

        def set_geometry(self, transverse_number_of_rows, longitudinal_number_of_rows,
                         transverse_pitch, longitudinal_pitch, L_ccrs, arrangement):
            """Set the geometric configuration of the tube bank.

            Args:
                transverse_number_of_rows (int): Number of tube rows transverse to flow.
                longitudinal_number_of_rows (int): Number of tube rows along flow direction.
                transverse_pitch (float): Transverse pitch between tubes in meters.
                longitudinal_pitch (float): Longitudinal pitch between tubes in meters.
                L_ccrs (float): Tube length at the cross-section.
                finned_tube_section (float): Length of the tube covered by fins (m).
                arrangement (arrangementType): Tube arrangement type.

            Notes:
                Derived properties such as diagonal pitch and total tube count
                are computed automatically.
            """
            self.geom_set = True

            # Row numbers
            self.transverse_number_of_rows      = transverse_number_of_rows
            self.longitudinal_number_of_rows    = longitudinal_number_of_rows
            self.total_number_of_tubes          = self.transverse_number_of_rows * self.longitudinal_number_of_rows

            # Pitches
            self.transverse_pitch               = transverse_pitch
            self.longitudinal_pitch             = longitudinal_pitch
            self.diagonal_pitch                 = np.sqrt((1/4)*self.transverse_pitch**2 + self.longitudinal_pitch**2)

            # Lengths
            self.L_ccrs                         = L_ccrs
            
            # Arrangement
            self.arrangement                    = arrangement 

        @cached_property
        def phi_parameter(self):
            """Compute the phi parameter used in flow area calculations.

            Returns:
                float: Ratio involving transverse and diagonal pitch minus
                the conventional diameter.
            """
            if not self.parent.Tube.geom_set or not self.parent.Fin.geom_set:
                raise ValueError("Tube and Fin geometries must be set before computing phi parameter.")
            
            d_cl = self.parent.Conventional_diameter
            phi_cl = (self.transverse_pitch-d_cl)/(self.diagonal_pitch-d_cl)
            return phi_cl

    def __init__(self):
        # Initialize inner classes
        self.Tube = self._Tube(self)
        self.Duct = self._Duct(self)
        self.Fin  = self._Fin(self)
        self.Bank = self._Bank(self)

    @cached_property
    def Conventional_diameter(self):
        """Compute the conventional finned-tube diameter.

        Returns:
            float: Conventional diameter based on fin geometry.
        """
        d_cl = self.Tube.outer_diameter + \
        ((2*self.Fin.fin_height*self.Fin.average_fin_thickness)/(self.Fin.fin_spacing))
        return d_cl

    @cached_property
    def Free_flow_area(self):
        """Compute the free-flow area available for air passage.

        Returns:
            float: Free-flow area in square meters.

        Notes:
            The expression depends on the phi parameter and conventional diameter.
        """
        d_cl = self.Conventional_diameter
        phi_cl = self.phi_parameter

        if phi_cl <= 2:
            F = self.Duct.a*self.Duct.b - self.Bank.transverse_number_of_rows*self.Bank.L_ccrs*d_cl
        else:
            F = (self.Duct.a*self.Duct.b - self.Bank.transverse_number_of_rows*self.Bank.L_ccrs*d_cl)*(2/phi_cl)
        return F

    @cached_property
    def Equivalent_diameter(self):
        """Compute the equivalent hydraulic diameter for the finned geometry.

        Returns:
            float: Equivalent diameter in meters.

        Notes:
            Adjustments are applied depending on the phi parameter.
        """
        phi_cl = self.phi_parameter
        
        d_eq = (2*(self.Fin.fin_spacing*(self.Bank.transverse_pitch-self.Tube.outer_diameter)
                - 2*self.Fin.fin_height*self.Fin.average_fin_thickness)) / \
               (2*self.Fin.fin_height + self.Fin.fin_spacing)
        
        if phi_cl > 2:
            d_eq = (2*d_eq)/phi_cl

        return d_eq

    @cached_property
    def A_total_over_F(self):
        """Compute the ratio of total heat-transfer area to free-flow area.

        Returns:
            float: Ratio A_total / F.
        """
        return (np.pi*(self.Tube.outer_diameter*self.Fin.fin_spacing
                + 2*self.Fin.fin_height*self.Fin.average_fin_thickness
                + 2*self.Fin.fin_height*(self.Fin.fin_height+self.Tube.outer_diameter))) / \
               (self.Bank.transverse_pitch*self.Fin.fin_spacing
                - (self.Tube.outer_diameter*self.Fin.fin_spacing
                + 2*self.Fin.fin_height*self.Fin.average_fin_thickness))

    @cached_property
    def S_1_over_S_2(self):
        """Compute the ratio of transverse pitch to longitudinal pitch.

        Returns:
            float: S1 / S2.
        """
        return self.Bank.transverse_pitch/self.Bank.longitudinal_pitch

    @cached_property
    def sigma_1(self):
        """Compute the transverse pitch-to-diameter ratio.

        Returns:
            float: Transverse pitch divided by tube outer diameter.
        """
        return self.Bank.transverse_pitch / self.Tube.outer_diameter
    
    @cached_property
    def sigma_2(self):
        """Compute the longitudinal pitch-to-diameter ratio.

        Returns:
            float: Longitudinal pitch divided by tube outer diameter.
        """
        return self.Bank.longitudinal_pitch / self.Tube.outer_diameter

    @cached_property
    def A_r(self):
        """Compute the total fin heat-transfer area.

        Returns:
            float: Fin surface area in square meters.

        Notes:
            The expression depends on fin type and tube bank configuration.
        """
        match self.Fin.fin_type:
            case finType.Rectangular:
                A_r = 2 * (self.Fin.square_height**2 - 0.785 * self.Tube.outer_diameter**2
                        + 2 * self.Fin.square_height*self.Fin.average_fin_thickness) \
                        * self.Bank.finned_tube_segment / self.Fin.fin_spacing \
                        * self.Bank.total_number_of_tubes
            case finType.Spiral | finType.Circular:
                A_r = np.pi/2 * (self.Fin.finning_diameter**2 - self.Tube.outer_diameter**2
                        + 2 * self.Fin.finning_diameter * self.Fin.average_fin_thickness) \
                        * self.Bank.finned_tube_segment / self.Fin.fin_spacing \
                        * self.Bank.total_number_of_tubes
        return A_r
        
    @cached_property
    def A_t(self):
        """Compute the tube heat-transfer area. Excluding the area covered by fins.
        
        Returns:
            float: Tube surface area in square meters.
        """
        L_covered_by_fins = self.Bank.finned_tube_segment * (self.Fin.average_fin_thickness/self.Fin.fin_spacing)
        L_bare_tube = (self.Bank.L_ccrs - L_covered_by_fins)* self.Bank.total_number_of_tubes
        
        A_t = np.pi * self.Tube.outer_diameter * L_bare_tube
        return A_t

    @cached_property
    def A(self):
        """Compute the total heat-transfer area.

        Returns:
            float: Sum of fin area and tube area.
        """
        return self.A_r + self.A_t

class ControlVolume:
    pass

class system:
    pass

class Node: # Fluid state at a specific position. 
    def __init__(self,fluid_type,Input_1,Input_2):

        self.fluid  = pf.Fluid(fluid_type)
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

    # ----- pressure drop ----- #

    def pressure_drop(self):
        '''Determines the type of pressure drop based on the flow type and phase change'''
        phase_in    = self.Node_in.fluid.phase.name
        phase_out   = self.Node_out.fluid.phase.name

        if phase_in in ("Liquid", "Gas") and phase_out in ("Liquid", "Gas"):
            dp = self.monophase_pressure_drop(self.Node_in.fluid,self.Node_out.fluid)

        elif phase_in == "TwoPhase" and phase_out == "TwoPhase":
            dp = self.two_phase_pressure_drop()

        elif phase_in != phase_out:
            # Phase change direction
            # Representative pressure drop for phase change segment
            if phase_in == "Gas" and phase_out == "TwoPhase":
                dp = self.monophase_pressure_drop(self.Node_in.fluid,self.Node_in.fluid)
 
            elif phase_in == "TwoPhase" and phase_out == "Liquid":
                dp = self.monophase_pressure_drop(self.Node_out.fluid,self.Node_out.fluid)
        else:
            raise ValueError("Invalid phase change direction for pressure drop calculation")

        return dp
    
    def monophase_pressure_drop(self,Fluid_in,Fluid_out):

        # Calculate flow velocities
        u_in       = self.m_dot / (Fluid_in.density * self.Geometry.Tube.inner_area)
        u_out      = self.m_dot / (Fluid_out.density * self.Geometry.Tube.inner_area)

        # Calculate average properties
        rho_avg    = (Fluid_in.density + Fluid_out.density) / 2
        u_avg      = self.m_dot / (rho_avg * self.Geometry.Tube.inner_area)
        mu_avg     = (Fluid_in.dynamic_viscosity + Fluid_out.dynamic_viscosity) / 2           
        Re_avg     = self.Reynolds_number(rho_avg,u_avg,mu_avg,self.Geometry.Tube.inner_diameter)

        # Acceleration pressure drop
        DP_A        = 0.5 * (Fluid_out.density * u_out**2 - Fluid_in.density * u_in**2)
        
        # Frictional pressure drop (Darcy-Weisbach equation)
        # Calculate friction factor depending on flow regime
        zeta_fr = self.friction_factor(Re_avg)
        
        # Calculate frictional pressure drop
        DP_R    = zeta_fr * (self.L/self.Geometry.Tube.inner_diameter) * 0.5 * (rho_avg * u_avg**2)
        
        # Calculate gravitational pressure drop
        DP_G    = 0 # Assuming horizontal flow for now

        # Total pressure drop of internal segment
        DP_tot  = DP_A + DP_R + DP_G

        return(DP_tot)
    
    def two_phase_pressure_drop(self):

        Node_in_liq     = self.Node_in.liquid_phase()
        Node_out_liq    = self.Node_out.liquid_phase()
        Node_in_vap     = self.Node_in.vapor_phase()
        Node_out_vap    = self.Node_out.vapor_phase()

        # Calculate average flow velocities
        u_liq_avg   = self.internal_flow_velocity(self.m_dot,(Node_in_liq.density + Node_out_liq.density)/2,self.Geometry.Tube.inner_area)
        u_vap_avg   = self.internal_flow_velocity(self.m_dot,(Node_in_vap.density + Node_out_vap.density)/2,self.Geometry.Tube.inner_area)

        # Calculate average properties
        rho_avg_liq = (Node_in_liq.density + Node_out_liq.density) / 2
        rho_avg_vap = (Node_in_vap.density + Node_out_vap.density) / 2
        mu_avg_liq  = (Node_in_liq.dynamic_viscosity + Node_out_liq.dynamic_viscosity) / 2
        mu_avg_vap  = (Node_in_vap.dynamic_viscosity + Node_out_vap.dynamic_viscosity) / 2

        # Calculate Reynolds numbers for liquid and vapor phases
        Re_avg_liq  = self.Reynolds_number(rho_avg_liq,u_liq_avg,mu_avg_liq,self.Geometry.Tube.inner_diameter)
        Re_avg_vap  = self.Reynolds_number(rho_avg_vap,u_vap_avg,mu_avg_vap,self.Geometry.Tube.inner_diameter)

        # average quality
        x_avg       = (self.Node_in.fluid.quality + self.Node_out.fluid.quality) / 2

        # Void fraction:
        eps_in      = self.smith_void_fraction(self.Node_in.fluid.quality,Node_in_vap.density,Node_in_liq.density)
        eps_out     = self.smith_void_fraction(self.Node_out.fluid.quality,Node_out_vap.density,Node_out_liq.density)

        # Accelerational pressure drop

        # Homogeneous model:
        P_A_in      = ((self.Node_in.fluid.quality**2)/(Node_in_vap.density * eps_in)+((1-self.Node_in.fluid.quality)**2)/(Node_in_liq.density*(1-eps_in)))
        P_A_out     = ((self.Node_out.fluid.quality**2)/(Node_out_vap.density * eps_out)+((1-self.Node_out.fluid.quality)**2)/(Node_out_liq.density*(1-eps_out)))
        DP_A        = (self.m_dot/self.Geometry.Tube.inner_area)**2 * (P_A_out - P_A_in)

        # Heterogeneous model:
        # u_g_in      = self.internal_flow_velocity(self.m_dot*self.Node_in.fluid.quality,Node_in_vap.density,self.Geometry.Tube.inner_area*eps_in)
        # u_g_out     = self.internal_flow_velocity(self.m_dot*self.Node_out.fluid.quality,Node_out_vap.density,self.Geometry.Tube.inner_area*eps_out)
        # u_l_in      = self.internal_flow_velocity(self.m_dot*(1-self.Node_in.fluid.quality),Node_in_liq.density,self.Geometry.Tube.inner_area*(1-eps_in))
        # u_l_out     = self.internal_flow_velocity(self.m_dot*(1-self.Node_out.fluid.quality),Node_out_liq.density,self.Geometry.Tube.inner_area*(1-eps_out))
        # P_A_in      = (u_g_in**2 * Node_in_vap.density * eps_in + u_l_in**2 * Node_in_liq.density * (1-eps_in))
        # P_A_out     = (u_g_out**2 * Node_out_vap.density * eps_out + u_l_out**2 * Node_out_liq.density * (1-eps_out))
        # DP_A        = P_A_out - P_A_in
        # print(DP_A)

        # Friction factor for liquid and vapor phases
        zeta_fr_liq = self.friction_factor(Re_avg_liq)
        zeta_fr_vap = self.friction_factor(Re_avg_vap)
        M_dot       = self.m_dot / self.Geometry.Tube.inner_area

        # Pressure drop of liquid and vapor phases
        DP_R_liq    = zeta_fr_liq * (M_dot**2)/(2 * rho_avg_liq * self.Geometry.Tube.inner_diameter)
        DP_R_vap    = zeta_fr_vap * (M_dot**2)/(2 * rho_avg_vap * self.Geometry.Tube.inner_diameter)

        # Assigning coefficients
        A       = DP_R_liq
        B       = DP_R_vap
        G       = A + 2*(B-A)*x_avg
        C       = 3

        # Frictional pressure drop
        DP_R    = (G * (1-x_avg)**(1/C) + B*x_avg**C) * self.L

        # Calculate gravitational pressure drop
        DP_G    = 0 # Assuming horizontal flow for now

        # Total pressure drop of internal segment
        DP_tot  = DP_A + DP_R + DP_G        

        return(DP_tot)

    def friction_factor(self,Re):
        if 3000 < Re < 100000:
            zeta = 0.3614/((Re**(1/4)))
        elif 2*10**4 < Re < 2*10**6:
            zeta = 0.3614/((Re**(1/4)))
        elif Re > 10**6:
            f = lambda z: 1/z + 0.8 - 2*np.log(Re*np.sqrt(z))
            zeta = sp.optimize.fsolve(f, x0=0.001)[0]
        return zeta

    def internal_flow_velocity(self,m_dot,density,inner_area):
        return m_dot / (density * inner_area)

    def smith_void_fraction(self,quality,rho_g,rho_l):
        K       = 0.4
        eps     = (1+ (rho_g/rho_l) * K * (1/quality - 1) + (rho_g/rho_l) * (1-K) * (1/quality - 1) * ((rho_l/rho_g + K *(1/quality - 1))/(1 + K * (1/quality - 1)))**(1/2) )**(-1)
        return eps

    def zivi_void_fraction(self,quality,rho_g,rho_l):
        eps     = (1+((1-quality)/quality)*(rho_g/rho_l)**(2/3))**(-1)
        return eps

    # ----- heat transfer ----- #

    def heat_transfer(self):

        phase_in    = self.Node_in.fluid.phase.name
        phase_out   = self.Node_out.fluid.phase.name        

        if phase_in in ("Liquid", "Gas") and phase_out in ("Liquid", "Gas"):
            h = self.monophase_heat_transfer(self.Node_in.fluid,self.Node_out.fluid)

        elif phase_in == "TwoPhase" and phase_out == "TwoPhase":
            h = self.twophase_heat_transfer(self.Node_in,self.Node_out)

        elif phase_in != phase_out:
            # Phase change direction
            # Representative pressure drop for phase change segment
            if phase_in == "Gas" and phase_out == "TwoPhase":
                h = self.monophase_heat_transfer(self.Node_in.fluid,self.Node_in.fluid)

            elif phase_in == "TwoPhase" and phase_out == "Liquid":
                h = self.monophase_heat_transfer(self.Node_in.fluid,self.Node_in.fluid)

            else:
                raise ValueError("Error, boiling")
        return h

    def monophase_heat_transfer(self,Fluid_in,Fluid_out):
        # Average values:
        rho_avg     = (Fluid_in.density + Fluid_out.density) / 2
        u_avg       = self.m_dot / (rho_avg * self.Geometry.Tube.inner_area)
        mu_avg      = (Fluid_in.dynamic_viscosity + Fluid_out.dynamic_viscosity) / 2           
        cp_avg      = (Fluid_in.specific_heat + Fluid_out.specific_heat) / 2
        k_avg       = (Fluid_in.conductivity + Fluid_out.conductivity) / 2

        # Dimensionless number calculations:
        Re_avg      = self.Reynolds_number(rho_avg,u_avg,mu_avg,self.Geometry.Tube.inner_diameter)
        Pr_avg      = self.Prandtl_number(cp_avg,mu_avg,k_avg)

        n_tem       = 0.25
        C_tem       = (mu_avg/mu_avg)**n_tem # until wall temperature is decided

        # Colburn equation: 
        if 10000 > Re_avg:
            print("Error, Laminar flow")

        elif 10000 < Re_avg < 1*10**6 and 0.7 < Pr_avg < 2: 
            h       = 0.023 * k_avg/self.Geometry.Tube.inner_diameter * Re_avg**(0.8)*Pr_avg**(0.4) * C_tem
        
        elif 4*10**3 < Re_avg < 5*10**6 and 0.1 < Pr_avg < 2000:
            zeta    = 1+900/Re_avg
            lmbda   = (1.82*np.log(Re_avg)-1.64)**(-2)  
            h       = k_avg/self.Geometry.Tube.inner_diameter * ((0.125*zeta*Re_avg*Pr_avg*C_tem)/(lmbda+4.5*zeta**(0.5)*(Pr_avg**(0.666-1))))

        elif Re_avg > 10000 and 0.7 < Pr_avg < 160: # skal lige eftertjekkes
            h       = 0.023 * k_avg/self.Geometry.Tube.inner_diameter * Re_avg**(0.8)*Pr_avg**(0.4)

        else: 
            print("Error, out of correlation range")

        return h

    def twophase_heat_transfer(self,Node_in,Node_out):
    
        # Average fluid quality
        x_avg           = (Node_in.fluid.quality+Node_out.fluid.quality) / 2
        Node_in_liq     = Node_in.liquid_phase()
        Node_out_liq    = Node_out.liquid_phase()

        # Calculate average liquid velocity
        u_liq_avg   = self.m_dot / (( (Node_in_liq.density + Node_out_liq.density) / 2) * self.Geometry.Tube.inner_area)

        # Calculate average liquid properties
        rho_avg_liq = (Node_in_liq.density + Node_out_liq.density) / 2
        mu_avg_liq  = (Node_in_liq.dynamic_viscosity + Node_out_liq.dynamic_viscosity) / 2
        cp_avg_liq  = (Node_in_liq.specific_heat + Node_out_liq.specific_heat) / 2
        k_avg_liq   = (Node_in_liq.conductivity + Node_out_liq.conductivity) / 2

        # Average Re and Pr values for liquid only
        Re_avg_liq  = self.Reynolds_number(rho_avg_liq,u_liq_avg,mu_avg_liq,self.Geometry.Tube.inner_diameter)        

        if Re_avg_liq > 3000: 

            Pr_avg_liq  = self.Prandtl_number(cp_avg_liq,mu_avg_liq,k_avg_liq)

            # Heat transfer of liquid only 
            h_L         = 0.023*Re_avg_liq**0.8 *Pr_avg_liq**(0.4) *k_avg_liq/self.Geometry.Tube.inner_diameter

            # Critical pressure
            p_cr        = self.Node_in.fluid.critical_pressure
            
            # Average actual pressure
            p_avg       = (self.Node_in.fluid.pressure + self.Node_out.fluid.pressure) / 2

            # Reduced pressure
            p_r         = p_avg / p_cr
            
            # Two phase heat transfer coefficient
            h_TP        = h_L * ((1-x_avg)**(0.8)+(3.8*x_avg**(0.76)*(1-x_avg)**(0.04))/(p_r**(0.38)))        

        else:
            print("Error, laminar flow")

        return h_TP

    # ----- General functions for internal flow ----- #

    def Reynolds_number(self,rho,u,mu,d):
        Re = (rho * u * d) / mu
        return Re
    
    def Prandtl_number(self,cp,mu,k):
        Pr = (cp*mu)/k
        return Pr

class External: # External correlations
    def __init__(self,mass_flow,rownumber,node_in,node_out,geometry,flue_in):
        self.z_2        = rownumber # Row number in the bank
        self.node_in    = node_in   # Node object representing the fluid state
        self.node_out   = node_out  # Node object representing the fluid state
        self.geometry   = geometry  # Geometry class
        self.mass_flow  = mass_flow # Mass flow
        self.flue_in    = flue_in   # Flue gas inlet parameters

    def heat_transfer(self,T_outer_wall):
        F = self.geometry.Free_flow_area
        d = self.geometry.outer_diameter
        
        avg_density = (self.node_in.fluid.density + self.node_out.fluid.density) / 2
        avg_Pr = (self.node_in.Prandtl_number + self.node_out.Prandtl_number) / 2
        u_g = self.mass_flow/(F*avg_density)
        
        avg_nu = (self.node_in.kinematic_viscosity + self.node_out.kinematic_viscosity) / 2
        avg_thermal_conductivity = (self.node_in.conductivity + self.node_out.conductivity) / 2

        avg_Re = (u_g * d) / avg_nu

        match self.geometry.Bank.arrangement:
            case arrangementType.Inline:
                X = 4 * (2 + self.geometry.Psi_r/7 - self.geometry.sigma_2) 
            case arrangementType.Staggered:
                X = self.geometry.sigma_1 / self.geometry.sigma_2  - 1.26 / self.geometry.Psi_r - 2
            case _: 
                NotImplementedError("Only Inline and Staggered tube arrangements are implemented for external heat transfer coefficient calculation.")

        n = 0.7 + 0.08 * np.tah(X) + 0.005 * self.geometry.Psi_r
        C_q = (1.36-np.tanh(X))*( (1.1)/(self.geometry.Psi_r+8) - 0.014)

        # Calculate C_z based on number of rows and arrangement
        if self.geometry.Bank.transverse_number_of_rows >= 8:
            C_z = 1
        elif self.geometry.Bank.arrangement == arrangementType.Inline or (self.geometry.sigma_1/self.geometry.sigma_2) >= 2:
            C_z = 3.5 * self.geometry.Bank.transverse_number_of_rows**(0.03) - 2.72
        elif self.geometry.Bank.arrangement == arrangementType.Staggered and (self.geometry.sigma_1/self.geometry.sigma_2) < 2:
            C_z = 3.15 * self.geometry.Bank.transverse_number_of_rows**(0.05) - 2.5
        else:
            raise ValueError("Error in calculating C_z for external heat transfer coefficient.")
        
        # Determine external heat transfer coefficient:
        h_c = 1.13*C_z*C_q*avg_thermal_conductivity/d * avg_Re**n * avg_Pr**(1/3)
        
        
        beta = np.sqrt(2*h_c/(self.geometry.Fin.average_fin_thickness*))
        
        
        # Determine the reduced heat transfer coefficient:
        psi_e = 1 - 0.016 * (self.geometry.finning_diameter/d-1)*(1+np.tanh(2*beta*l_r-1))

    def pressure_drop(self):

        # Average properties
        rho_avg = (self.flue_in.fluid.density + self.node_out.fluid.density) / 2
        nu_avg  = (self.flue_in.fluid.kinematic_viscosity + self.node_out.fluid.kinematic_viscosity) / 2

        # Free flow area
        F       = self.geometry.Free_flow_area

        # Equivalent diameter
        d_eq    = self.geometry.Equivalent_diameter
        
        # Flue gas velocity
        u_g     = self.mass_flow/(F*rho_avg)

        # Flue gas Reynold number
        Re_eq   = (u_g*d_eq)/nu_avg
        
        if 5*10**3< Re_eq < 6 * 10**4:
            #Correction factor for small row numbers:
            if self.geometry.Bank.arrangement == 'staggered' and self.z_2 < 6:
                C_zm = np.exp(0.1*(6/self.z_2-1))

            elif self.geometry.Bank.arrangement == 'inline' and self.z_2 < 6:
                C_zm = 1+(0.65/(self.z_2**3))

            else:
                C_zm = 1

            # coefficient in similarity equation for aerodynamic resistance
            if self.geometry.Bank.arrangement == 'staggered':
                n = 0.17*(self.geometry.A_total_over_F)**(0.25) * (self.geometry.S_1_over_S_2)**(0.57) * np.exp(-0.36*(self.geometry.S_1_over_S_2))
                C_r = 2.8 * (self.geometry.A_total_over_F)**(0.53) * (self.geometry.S_1_over_S_2)**(1.30) * np.exp(-0.90*(self.geometry.S_1_over_S_2))
            elif self.geometry.Bank.arrangement == 'inline':
                if self.geometry.S_1_over_S_2 <= 2.1:
                    n = (self.geometry.A_total_over_F)**(0.08) * (0.184-0.088*(self.geometry.S_1_over_S_2))
                    C_r = 2.5 * (self.geometry.A_total_over_F)**(0.25) * np.exp(-1.70*(self.geometry.S_1_over_S_2))
                elif self.geometry.S_1_over_S_2 > 2.1:
                    n = 0
                    C_r = (self.geometry.A_total_over_F)**(0.10) * (0.132-0.016*(self.geometry.S_1_over_S_2))
            
            # Resistance coefficient is calculated:
            zeta_0 = C_zm * C_r * Re_eq**(-n)
    
            # Constant constant
            C_op = 1.1

            # Delta P calculation:
            delta_P = C_op * zeta_0 * self.z_2 * (rho_avg*u_g**2)/2 

        else:
            raise ValueError("Error, external flow out of correlation range.")

        return(delta_P)
        # Phi parameter:

    def pressure_drop_2(self):
        pass
    
    
class Control_Volume:
    def __init__(self,num,Internal,External):
        
        self.num = num 
        self.Internal = Internal
        self.External = External
        
    def solve_wall_temperature(self):
        T_oo_inner = (self.Internal.Node_in.fluid.temperature + self.Internal.Node_out.fluid.temperature) / 2
        T_oo_external = (self.External.node_in.fluid.temperature + self.External.node_out.fluid.temperature) / 2
        
        




if __name__ == "__main__":
    # Geom definition
    Geom = Geometry()
    Geom.Tube.set_geometry(0.028,0)
    Geom.Duct.set_geometry(0.56,0.5)
    Geom.Fin.set_geometry(8e-4,0.0135,0.003)
    Geom.Bank.set_geometry(9,30,0.059,0.051,0.5,'staggered')
    
    P_in_air = (-4e3) + 966e2 
    T_in_air = 230
    
    P_in_water = 135e5 # bar to Pa
    T_in_water = 360
    
    # phase change check
    Node_1 = Node(FluidsList.Water,1,pf.Input.temperature(340),pf.Input.pressure(13500000))
    Node_2 = Node(FluidsList.Water,1,pf.Input.quality(0.9),pf.Input.pressure(13500000))
    
    Node_1 = Node(FluidsList.Water,1,pf.Input.temperature(340),pf.Input.pressure(13500000))
    Node_2 = Node(FluidsList.Water,1,pf.Input.quality(1),pf.Input.pressure(13500000))
    
    Node_internal_S1 = Node(FluidsList.Water,0,pf.Input.temperature(T_in_water),pf.Input.pressure(P_in_water))
    Node_internal_S2 = Node(FluidsList.Water,0,pf.Input.temperature(T_in_water-10),pf.Input.pressure(P_in_water))
    
    Node_external_S1 = Node(FluidsList.Air,1,pf.Input.temperature(T_in_air),pf.Input.pressure(P_in_air))
    Node_external_S2 = Node(FluidsList.Air,1,pf.Input.temperature(T_in_air+10),pf.Input.pressure(P_in_air))
    
    Internal_1 = Internal(0.5,Node_internal_S1,Node_internal_S2,Geom,0.5)
    print(Internal_1.heat_transfer())
    Internal_2 = Internal(0.5,Node_1,Node_2,Geom,0.5)
    print(Internal_2.heat_transfer())
    External = External(0.5,1,Node_external_S1,Node_external_S2,Geom)
    
    
    
    
    


