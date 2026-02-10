import numpy as np
import scipy as sp
import re
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from scipy.optimize import fsolve, minimize, root
import plotly.graph_objects as go
from plotly.colors import sample_colorscale

class node: 
    def __init__(self, num, temperature,type):
        self.num = num
        self.temperature = temperature
        self.type = type  # 'hs' or 'fg'
    def __repr__(self):
        return f"Node {self.num} of type {self.type}"


class control_volume:
    def __init__(self,num):
        self.num = num
        self.is_boundary_fg = False
        self.is_boundary_hs = False
        
        # Calculation coefficients
        N_cv = N * M  # Number of control volumes 

        # Hard coded for simplicity
        self.UA = 40 / N_cv 
        self.m_dot_hs = 0.05 / N
        self.m_dot_fg = 0.5 / M 

        self.c_p_hs = 1005 # J/kgK
        self.c_p_fg = 100 # J/kgK

    def set_nodes(self, node_hs_in, node_hs_out, node_fg_in, node_fg_out):
        self.node_hs_in = node_hs_in
        self.node_hs_out = node_hs_out
        self.node_fg_in = node_fg_in
        self.node_fg_out = node_fg_out
        
    def set_boundary_hs(self,T_hs_BC):
        self.is_boundary_hs = True
        self.node_hs_in.temperature = T_hs_BC
        
    def set_boundary_fg(self,T_fg_BC):
        self.is_boundary_fg = True
        self.node_fg_in.temperature = T_fg_BC
    
    def calculate(self):
        # Energy balance for hot side
        Q_hs = self.m_dot_hs * self.c_p_hs * (self.node_hs_in.temperature - self.node_hs_out.temperature)
        # Energy balance for flue gas side
        Q_fg = self.m_dot_fg * self.c_p_fg * (self.node_fg_out.temperature - self.node_fg_in.temperature)
        # Heat exchanged
        Q_out = self.UA * ((self.node_hs_in.temperature + self.node_hs_out.temperature)/2 - (self.node_fg_in.temperature + self.node_fg_out.temperature)/2)
        return Q_hs, Q_fg, Q_out    
    
    # Properties to get min and max temperatures for visualization
    @property
    def min_temp_hs(self):
        return min(self.node_hs_in.temperature, self.node_hs_out.temperature)
    @property
    def max_temp_hs(self):
        return max(self.node_hs_in.temperature, self.node_hs_out.temperature)
    @property
    def min_temp_fg(self):
        return min(self.node_fg_in.temperature, self.node_fg_out.temperature)
    @property
    def max_temp_fg(self):
        return max(self.node_fg_in.temperature, self.node_fg_out.temperature)
    
    def __repr__(self):
        return f"Control Volume {self.num}"
    
class system:
    def __init__(self, M, N):
        self.M = M
        self.N = N
        
        # Create control volumes
        # Control volumes are numbered column wise 
        # CV[M][N] where M is number of control volumes in a row and N is number of rows
            # So for a 2x2 system control volume 0 is at (0,0), control volume 1 is at (1,0), control volume 2 is at (0,1) etc.  
        self.control_volumes = np.array([[control_volume(i*M+j) for j in range(M)] for i in range(N)])
        
        # Create nodes and assign to control volumes
        self._create_nodes()
    
    def _create_nodes(self):
        '''
        Create nodes for the system and assign them to control volumes. Each output node is the input node for the next control volume.
        For the flue gas side, the same node is shared between control volumes in the same row.
        '''
        
        hs_node_num = 0
        fg_node_num = 0
        
        for row_control_volume in self.control_volumes:
            if row_control_volume[0].num == 0:
                # First control volume in the system create both inlet and outlet nodes
                node_fg_in = node(fg_node_num, None, 'fg')
                fg_node_num += 1
                node_fg_out = node(fg_node_num, None, 'fg')
                fg_node_num += 1
            
            for control_volume in row_control_volume:
                if control_volume.num == 0:
                    # First control volume in the system create both inlet and outlet nodes
                    node_hs_in = node(hs_node_num, None, 'hs')
                    hs_node_num += 1
                    node_hs_out = node(hs_node_num, None, 'hs')
                    hs_node_num += 1
                
                # set nodes
                control_volume.set_nodes(node_hs_in, node_hs_out, node_fg_in, node_fg_out)
                
                # Prepare nodes for next control volume
                node_hs_in = node_hs_out
                node_hs_out = node(hs_node_num, None, 'hs')
                hs_node_num += 1
            
            # Prepare nodes for next row of control volumes
            node_fg_in = node_fg_out
            node_fg_out = node(fg_node_num, None, 'fg')
            fg_node_num += 1
            
    
    def set_boundary_conditions(self, T_hs_BC, T_fg_BC):
        '''
        Set boundary conditions for the system. T_hs_BC is the inlet temperature for the hot side. T_fg_BC is the inlet temperature for the flue gas side.
        '''
        # Set boundary conditions for hot side
        self.control_volumes[0,0].set_boundary_hs(T_hs_BC)
        # Set boundary conditions for flue gas side
        for control_volume in self.control_volumes[0,:]:
            control_volume.set_boundary_fg(T_fg_BC)
    
    def update_temp(self,T_guess):
        '''
        Updates the temperatures of the nodes in the system based on the provided temperature guesses.
        T_guess is a concatenated list of T_hs and T_fg guesses without boundary conditions.
        '''
        # Unpack temperature guesses
        
        # There is always the same number of hot side temperatures as control volumes
        # Since for the first row the inlet temperature is known from BCs
        # And for subsequent rows the outlet temperature of the previous row is the inlet temperature
        T_hs_guess = T_guess[:(self.M*self.N)]
        T_fg_guess = T_guess[(self.M*self.N):]
        
        # Update hot side temperatures
        for control_volume in self.control_volumes.flatten():
            if control_volume.is_boundary_hs:
                control_volume.node_hs_out.temperature = T_hs_guess[control_volume.num]
            else:
                control_volume.node_hs_in.temperature = T_hs_guess[control_volume.num-1]
                control_volume.node_hs_out.temperature = T_hs_guess[control_volume.num]
            
        # Update flue gas side temperatures
        for row_num,row in enumerate(self.control_volumes):
            for control_volume in row:
                if control_volume.is_boundary_fg or row_num == 0:
                    control_volume.node_fg_out.temperature = T_fg_guess[row_num]
                else:
                    control_volume.node_fg_in.temperature = T_fg_guess[row_num-1]
                    control_volume.node_fg_out.temperature = T_fg_guess[row_num]

    def solve_system(self,T_hs_BC,T_fg_BC):
        '''
        Solve the system of equations to find the temperatures of the nodes.
        '''
        
        self.set_boundary_conditions(T_hs_BC, T_fg_BC)
        # Initial guess for temperatures
        T_hs_guess = np.array([T_hs_BC - i*4 for i in range(self.M*self.N)])
        T_fg_guess = np.array([T_fg_BC + i*5 for i in range(self.N)])
        
        T_guess = np.concatenate((T_hs_guess, T_fg_guess))
        
        def equations(T_guess):
            self.update_temp(T_guess)
            eqs_hs = []
            eqs_fg = []
            for row in self.control_volumes:
                fg_list = []
                for control_volume in row:
                    Q_hs, Q_fg, Q_out = control_volume.calculate()
                    eqs_hs.append(Q_hs - Q_out)
                    fg_list.append(Q_fg - Q_out)
                
                eqs_fg.append(sum(fg_list))
            return np.array(eqs_hs + eqs_fg)
        
        sol = root(equations, T_guess, method='hybr')
        return sol
    
    # Functions to display the system
    def show_system(self):
        # Create a copy to manipulate for display
        control_volumes = self.control_volumes.copy()
        
        # Flip every second column to represent snaking flow path
        control_volumes[1::2, :] = np.flip(control_volumes[1::2, :], axis=1)
                
        print("System Control Volumes:")
        for row in control_volumes:
            print(*row,sep='\t')
            
    def show_nodes(self):
        '''
        Show all nodes in the system. in an array format with snaking flow path.
        '''
        # Create a copy to manipulate for display
        control_volumes = self.control_volumes.copy()
        # Flip every second column to represent snaking flow path
        control_volumes[1::2, :] = np.flip(control_volumes[1::2, :], axis=1)
        
        print("\nSystem Nodes:")
        for row_idx, row in enumerate(control_volumes):
            row_matrix = []
            swap_hs = (row_idx % 2 == 1)
            for control_volume in row:
                hs_left = control_volume.node_hs_in.num
                hs_right = control_volume.node_hs_out.num
                if swap_hs:
                    hs_left, hs_right = hs_right, hs_left
                bold_cv = f"\x1b[4m{control_volume.num}\x1b[0m"
                
                # Create 3x3 matrix of nodes for each control volume to represent the connections
                # Layout: fg_in (top), hs_in (left), CV (center), hs_out (right), fg_out (bottom)
                matrix = np.array([["", control_volume.node_fg_in.num, ""],
                                  [hs_left, bold_cv, hs_right],
                                  ["", control_volume.node_fg_out.num, ""]])
                
                # Concatenate matrices horizontally
                row_matrix.append(matrix)
            for i in range(len(row_matrix[0])):
                line_cells = [cell for block in row_matrix for cell in block[i]]
                string_cells = ["" if cell == "" else str(cell) for cell in line_cells]
                cell_lengths = [len(re.sub(r"\x1b\[[0-9;]*m", "", cell)) for cell in string_cells]
                max_width = max(cell_lengths, default=1)
                line = " ".join(cell.rjust(max_width) if cell else "".rjust(max_width) for cell in string_cells)
                print(line.rstrip())
                
    def show_temperatures(self):
        '''
        Show all the temperatures of the nodes in the system. in an array format with snaking flow path.
        '''
        # Create a copy to manipulate for display
        control_volumes = self.control_volumes.copy()
        # Flip every second column to represent snaking flow path
        control_volumes[1::2, :] = np.flip(control_volumes[1::2, :], axis=1)
        
        print("\nSystem Node Temperatures:")
        for row_idx, row in enumerate(control_volumes):
            row_matrix = []
            swap_hs = (row_idx % 2 == 1)
            for control_volume in row:
                hs_left_temp = control_volume.node_hs_in.temperature
                hs_right_temp = control_volume.node_hs_out.temperature
                if swap_hs:
                    hs_left_temp, hs_right_temp = hs_right_temp, hs_left_temp
                bold_cv = f"\x1b[4m{control_volume.num}\x1b[0m"
                
                # Create 3x3 matrix of temperatures for each control volume to represent the connections
                # Layout: fg_in (top), hs_in (left), CV (center), hs_out (right), fg_out (bottom)
                matrix = np.array([["", f"{control_volume.node_fg_in.temperature:.1f}", ""],
                                  [f"{hs_left_temp:.1f}", bold_cv, f"{hs_right_temp:.1f}"],
                                  ["", f"{control_volume.node_fg_out.temperature:.1f}", ""]])
                
                # Concatenate matrices horizontally
                row_matrix.append(matrix)
            for i in range(len(row_matrix[0])):
                line_cells = [cell for block in row_matrix for cell in block[i]]
                string_cells = ["" if cell == "" else str(cell) for cell in line_cells]
                cell_lengths = [len(re.sub(r"\x1b\[[0-9;]*m", "", cell)) for cell in string_cells]
                max_width = max(cell_lengths, default=1)
                line = " ".join(cell.rjust(max_width) if cell else "".rjust(max_width) for cell in string_cells)
                print(line.rstrip())

    def visualize_system(self, show_temps=False):
        """
        Shows a visual representation of the system using matplotlib. 
        Each control volume is represented as a rectangle with a vertical gradient representing the flue gas temperature and a colored line representing the hot side temperature. 
        The nodes are plotted as points with colors corresponding to their temperatures. The flow path is represented by the arrangement of the control volumes and lines.
        
        The nodes can be clicked to show their properties in an annotation. 
        """
        import matplotlib.pyplot as plt
        import matplotlib.colors as mcol
        import matplotlib.collections as mcollections
                
        def colored_line(x, y, c, ax,stroke=False,stroke_color='black', **lc_kwargs):
            """
            Plot a line with a color specified along the line by a third value.
            Parameters
            ----------
            x,y : array-like 
                The x and y coordinates of the line points.
            c : array-like
                The values used to color the line. Must be the same length as x and y.
            ax : matplotlib axis
                The axis to plot on.
            stroke : bool, optional
                Whether to add a stroke (outline) to the line for better visibility. Default is False.
            stroke_color : str, optional
                The color of the stroke if stroke is True. Default is 'black'.
            **lc_kwargs : keyword arguments
                Additional keyword arguments passed to LineCollection.
            """ 
            # Based on https://matplotlib.org/stable/gallery/lines_bars_and_markers/multicolored_line.html with modifications for stroke
            
            # Default the capstyle to butt so that the line segments smoothly line up
            default_kwargs = {"capstyle": "butt"}
            default_kwargs.update(lc_kwargs)
            
            # Compute the midpoints of the line segments. Include the first and last points
            # twice so we don't need any special syntax later to handle them.
            x = np.asarray(x)
            y = np.asarray(y)
            x_midpts = np.hstack((x[0], 0.5 * (x[1:] + x[:-1]), x[-1]))
            y_midpts = np.hstack((y[0], 0.5 * (y[1:] + y[:-1]), y[-1]))
            
            # Create a line below the main line to create a stroke effect
            if not stroke == False:
                ax.plot(x, y, color=stroke_color, linewidth=lc_kwargs.get('linewidth', 2) + stroke, zorder=1)        
            
            # Create line segments for the LineCollection. The segments are defined by their start, midpoint, and end coordinates.
            coord_start = np.column_stack((x_midpts[:-1], y_midpts[:-1]))[:, np.newaxis, :]
            coord_mid = np.column_stack((x, y))[:, np.newaxis, :]
            coord_end = np.column_stack((x_midpts[1:], y_midpts[1:]))[:, np.newaxis, :]
            segments = np.concatenate((coord_start, coord_mid, coord_end), axis=1)

            lc = mcollections.LineCollection(segments, **default_kwargs)
            lc.set_array(c)  # set the colors of each segment

            return ax.add_collection(lc)
        
        def vertical_gradient_rect(control_volume, horizontal_span, vertical_span, ax, cmap, norm):
            """
            Create a vertical gradient rectangle at the specified x and y coordinates with the given color values, colormap, and normalization.
            
            Parameters
            ----------
            control_volume : control_volume
                The control volume to plot.
            horizontal_span : float
                The horizontal span of the rectangle.
            vertical_span : float
                The vertical span of the rectangle.
            ax : matplotlib axis
                The axis to plot on.
            cmap : matplotlib colormap
                The colormap to use for coloring the rectangle.
            norm : matplotlib Normalize
                The normalization to apply to the color values.
            ax : matplotlib axis
                The axis to plot on.
            cmap : matplotlib colormap
                The colormap to use for coloring the rectangle.
            norm : matplotlib Normalize
                The normalization to apply to the color values.
            """
            x, y = control_volume.center
            
            # Flue gas gradient
            x = [x - horizontal_span/2, x + horizontal_span/2]
            y = [y - vertical_span/2, y + vertical_span/2]
            
            c = [control_volume.node_fg_in.temperature, control_volume.node_fg_out.temperature]
            
            return ax.imshow(
                [[c[0], c[0]], [c[1], c[1]]], # The color values for the rectangle vertices. 
                cmap=cmap,
                norm=norm,
                interpolation='gaussian', # Smooth and fast interpolation
                extent=[x[0], x[1], y[0], y[1]]  # The extent of the rectangle
            )
            
        def plot_nodes(control_volume, horizontal_span,vertical_span, ax, cm_hs, norm_hs, cm_fg, norm_fg, show_temps=False,is_first_last=False):
            '''
            Plot the nodes of a control volume on the given axis.
            
            Parameters
            ----------
            control_volume : control_volume
                The control volume to plot.
            horizontal_span : float
                The horizontal span of the nodes.
            vertical_span : float
                The vertical span of the nodes.
            ax : matplotlib axis
                The axis to plot on.
            cm_hs : matplotlib colormap
                The colormap to use for coloring the hot side nodes.
            norm_hs : matplotlib Normalize
                The normalization to apply to the color values for the hot side nodes.
            cm_fg : matplotlib colormap
                The colormap to use for coloring the flue gas nodes.
            norm_fg : matplotlib Normalize
                The normalization to apply to the color values for the flue gas nodes.
            show_temps : bool, optional
                Whether to display the temperatures on the nodes. Default is False.
            is_first_last : bool, optional
                Whether the control volume is the first or last in the row. This is used to determine whether to plot the inlet node. Default is False. 1 is first and -1 is last.
            '''
            
            x,y = control_volume.center
            
            def plot_node(control_volume, node, x, y, ax, cm, norm, show_temp=False):
                
                color = cm(norm(node.temperature))
                
                artist = ax.scatter(x,y, color=color, zorder=5,picker=True)  # Hot side outlet node
                artist.obj = node  # Attach the node object to the artist for later retrieval
                artist.obj.owner_num = control_volume.num
            
            # Plot hot side nodes
            if not is_first_last == False: 
                # Plot both inlet and outlet nodes
                plot_node(control_volume, control_volume.node_hs_in, x - horizontal_span/2, y, ax, cm_hs, norm_hs, show_temps)
                plot_node(control_volume, control_volume.node_hs_out, x + horizontal_span/2, y, ax, cm_hs, norm_hs, show_temps)
            else:
                # Plot only the outlet node since the inlet node is the same as the previous control volume's outlet node
                plot_node(control_volume, control_volume.node_hs_out, x + horizontal_span/2, y, ax, cm_hs, norm_hs, show_temps)
        
            # Plot flue gas nodes
            plot_node(control_volume, control_volume.node_fg_out, x, y - vertical_span/2, ax, cm_fg, norm_fg, show_temps)  
            if control_volume.is_boundary_fg:
                plot_node(control_volume, control_volume.node_fg_in, x, y + vertical_span/2, ax, cm_fg, norm_fg, show_temps)
        
        
        def plot_line_end(control_volume,horizontal_span,vertical_span,plot_length,direction,ax):
            '''
            Plot the line connecting the hot side inlet and outlet nodes of a row.
            
            Parameters
            ----------
            control_volume : control_volume
                The control volume to plot.
            horizontal_span : float
                The horizontal span of the control volume.
            vertical_span : float
                The vertical span of the control volume.
            plot_length : float
                The length of the line to plot.
            direction : int
                The direction of the line end. 1 for right, -1 for left.
            ax : matplotlib axis
                The axis to plot on.
            '''
            x, y = control_volume.center
            
            point1 = (x + direction*horizontal_span/2, y)
            point2 = (x + direction*horizontal_span/2 + direction*plot_length, y)
            point3 = (x + direction*horizontal_span/2 + direction*plot_length, y - vertical_span)
            point4 = (x + direction*horizontal_span/2 , y - vertical_span)
            
            ax.plot([point1[0], point2[0]], [point1[1], point2[1]], color='black', linewidth=2, zorder=4,linestyle='--')
            ax.plot([point2[0], point3[0]], [point2[1], point3[1]], color='black', linewidth=2, zorder=4,linestyle='--')
            ax.plot([point3[0], point4[0]], [point3[1], point4[1]], color='black', linewidth=2, zorder=4,linestyle='--')
            

                
        def plot_control_volume(control_volume, x, y, vertical_span, horizontal_span, ax, cm_hs, cm_fg, norm_hs, norm_fg,draw_boundary=True,is_first_last=False, show_temps=False):
            '''
            Plot a single control volume on the given axis with the specified colormaps and normalizations for the hot side and flue gas side temperatures.
            
            Parameters
            ----------
            control_volume : control_volume
                The control volume to plot.
            x : scalar
                The x coordinate of the control volume center.
            y : scalar
                The y coordinate of the control volume center.
            vertical_span : scalar
                The vertical span of the control volume.
            horizontal_span : scalar
                The horizontal span of the control volume.
            ax : matplotlib axis
                The axis to plot on.
            cm_hs : matplotlib colormap
                The colormap to use for the hot side temperatures.
            cm_fg : matplotlib colormap
                The colormap to use for the flue gas side temperatures.
            norm_hs : matplotlib Normalize
                The normalization to apply to the hot side temperature values.
            norm_fg : matplotlib Normalize
                The normalization to apply to the flue gas side temperature values.
            draw_boundary : bool, optional
                Whether to draw the boundary of the control volume. Default is True.
            is_first_last : bool, optional
                Whether the control volume is the first or last in the sequence. Default is False. 1 is first and -1 is last.
            show_temps : bool, optional
                Whether to display the temperatures on the control volume. Default is False.
            '''
            center = (x, y)
            
            control_volume.center = center # Store the center coordinates in the control volume for later use

            # Draw the flue gas side gradient rectangle first so that the hot side line and nodes are on top of it
            vertical_gradient_rect(control_volume, horizontal_span, vertical_span, ax, cmap=cm_fg, norm=norm_fg)
            
            # Draw a dashed boundary around the control volume for visibility
            if draw_boundary:
                bottom_left_corner = (x - horizontal_span/2, y - vertical_span/2)
                rect = patches.Rectangle(bottom_left_corner, horizontal_span, vertical_span, linewidth=1, edgecolor='black', facecolor='none', zorder=3,linestyle='--')
                ax.add_patch(rect)
            
            # Hot side colored line
            # Plot the hot side line with colors corresponding to the temperatures
            x_hs_min = x - horizontal_span/2
            x_hs_max = x + horizontal_span/2
            
            M = self.M
            Num_x_values = 100 / M # Number of x values for the hot side line
            
            x_hs = np.linspace(x_hs_min, x_hs_max, int(Num_x_values))
            y_hs = np.ones_like(x_hs) * y
            
            c_hs = np.linspace(control_volume.node_hs_in.temperature, control_volume.node_hs_out.temperature, len(x_hs))
            
            # Plot the hot side line with colors corresponding to the temperatures
            colored_line(x_hs, y_hs, c_hs, ax, stroke=3, stroke_color='black', linewidth=4, cmap=cm_hs, norm=norm_hs)
            # Plot the nodes on top of the hot side line so they are visible
            plot_nodes(control_volume, length_cv, height_cv, ax, cm_hs, norm_hs, cm_fg, norm_fg, show_temps=show_temps,is_first_last=is_first_last)
        
        # Create a copy to manipulate for display
        control_volumes = self.control_volumes.copy()
        
        # Flip every second column to represent snaking flow path
        control_volumes[1::2, :] = np.flip(control_volumes[1::2, :], axis=1)
        
        # Find max and min hot side and flue gas side temperatures for colormap normalization
        min_temp_hs = min(cv.min_temp_hs for cv in control_volumes.flatten())
        max_temp_hs = max(cv.max_temp_hs for cv in control_volumes.flatten())
        min_temp_fg = min(cv.min_temp_fg for cv in control_volumes.flatten())
        max_temp_fg = max(cv.max_temp_fg for cv in control_volumes.flatten())
        
        # Create colormaps and normalizations
        cm_hs = mcol.LinearSegmentedColormap.from_list("MyCmapName",["tab:blue","tab:red"])
        cm_fg = mcol.LinearSegmentedColormap.from_list("MyCmapName",["xkcd:grey blue","xkcd:grey pink"])
        
        norm_hs = mcol.Normalize(vmin=min_temp_hs, vmax=max_temp_hs)
        norm_fg = mcol.Normalize(vmin=min_temp_fg, vmax=max_temp_fg)
        
        # Create plot objects
        fig,ax = plt.subplots()
        
        # Number of control volumes in a row
        M = self.M
        
        # Number of longitudinal rows
        N = self.N
        
        plot_length = 10 # length of a tube row in meters should pull from geometry class 
        plot_height = 5
        
        length_cv = plot_length / M
        height_cv = plot_height / N
        
        line_end_length = plot_length / 20
        
        # Plot control volumes
        y = height_cv / 2 + height_cv * N  # y coordinates for the center of each control volume row
        final_cv_num = self.control_volumes.flatten()[-1].num
        
        for row in control_volumes:
            x = length_cv / 2
            # Determine flow direction for hot side
            if row[0].num < row[-1].num: 
                direction = 1
            else:
                direction = -1
                            
            for indx, control_volume in enumerate(row):
                # If first or last control volume in the row, 
                if indx == 0:
                    plot_control_volume(control_volume, x, y, height_cv, length_cv, ax, cm_hs, cm_fg, norm_hs, norm_fg, draw_boundary=True, is_first_last=1, show_temps=show_temps)
                elif indx == len(row) - 1:
                    plot_control_volume(control_volume, x, y, height_cv, length_cv, ax, cm_hs, cm_fg, norm_hs, norm_fg, draw_boundary=True, is_first_last=-1, show_temps=show_temps)
                else:
                    plot_control_volume(control_volume, x, y, height_cv, length_cv, ax, cm_hs, cm_fg, norm_hs, norm_fg, draw_boundary=True, show_temps=show_temps)
                    
                # Draw line ends for the last control volumes in the row if it's not the final control volume in the system
                if indx == 0 and direction == -1 and not control_volume.num == final_cv_num:
                    plot_line_end(control_volume,length_cv,height_cv,line_end_length,direction,ax)
                elif indx == len(row) - 1 and direction == 1 and not control_volume.num == final_cv_num:
                    plot_line_end(control_volume,length_cv,height_cv,line_end_length,direction,ax)
                
                
                x = x + length_cv
            y = y - height_cv

        # Define the onclick event handler for picking nodes
        Current_annotation = None  # Variable to keep track of the current annotation
        def onclick(event):
            print(f"artist picked: {event.artist}")
            print(f"xy {event.artist.get_offsets()}")  # Get the coordinates of the picked point
            print(f"node object: {event.artist.obj.temperature}")  # Access the attached node object and print it
            
            # Remove the previous annotation if it exists
            nonlocal Current_annotation
            if Current_annotation is not None:
                Current_annotation.remove() 
                Current_annotation = None
            
            # Create a new annotation for the picked node
            an = ax.annotate(
                f"Temp: {event.artist.obj.temperature:.1f} C \nNode: {event.artist.obj.num}, Owner CV: {event.artist.obj.owner_num} \n Type: {event.artist.obj.type}", 
                xy=(event.artist.get_offsets()[0]), xytext=(10,10),
                textcoords='offset points', 
                bbox=dict(boxstyle="round", fc="w"), 
                arrowprops=dict(arrowstyle="->"),
                zorder=10)
            
            # Store the current annotation so it can be removed when the next node is picked
            Current_annotation = an 
            # Redraw the canvas to show the annotation immediately
            fig.canvas.draw()
            
        
        # Connect the onclick event handler to the figure
        cid = fig.canvas.mpl_connect('pick_event', onclick)
        # Display the plot  
        plt.show()
        
            
if __name__ == "__main__":
    
    # Number of control volumes in a row
    M = 50
    # Number of longitudinal rows
    N = 4

    T_hs_BC = 400  # Inlet temperature hot side
    T_fg_BC = 200  # Inlet temperature flue gas side
    
    sys = system(M,N)
    sys.set_boundary_conditions(T_hs_BC, T_fg_BC)

    sol = sys.solve_system(T_hs_BC,T_fg_BC)
    sys.update_temp(sol.x)
    sys.visualize_system(show_temps=True)
    
    
    