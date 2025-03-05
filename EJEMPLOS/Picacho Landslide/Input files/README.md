# .DAT input file:

It is an input file where the necessary information to model the problem is stored.

**nline**: indicate the lines of text describing the problem. (1 by default)

**SWalg**: The algorithm used for the Shallow Water computations in the one representing landslides (Type of SW Algorithm). ic_SWAlg algorithm for SW computations. There are several options (0 by default, Professional users can use other types):

**nhist**: Number of histograms (0 by default, Professional users can only use this feature)

**ndimn**: dimension of the problem 1,2 3, Problem dimensions: A depth integrated two dimensional problem. (2 by default, Professional users can use other dimensions)

**soil**: Indicates the presence of soil particles, takes value 1 if affirmative and 0 if it does not proceed. (1 by default)

**water**: Indicates the presence of water, takes value 1 if affirmative and 0 if it does not proceed. (0 by default, Professional users can only use this feature)

**vps**: Indicates that virtual particles are considered in the analysis, it takes value 1 if affirmative and 0 if it does not proceed.(0 by default, Professional users can only use this feature)

**abss**: Indicates the activation of boundary conditions. (0 by default, Professional users can only use this feature)

**icunk**: It indicates the way that the file.pts has to be read indicate the way the landslide is described (indicate the way the landslide is described) and a limit for landslide depth. (6 by default, Professional users can use other types)

**h_inf_SW**: Unstable material's minimum height that can be considered by the code

**file.pts**: The name of .pts file

**pa_sph**: SPH algorithm for particle approximation (pa_sph) (2 by default, Professional users can use other types)

**nnps**: Nearest neighbor particle searching (nnps) method (2 by default, Professional users can use other types)

**sle**: smoothing length evolution (sle) algorithm (1 by default, Professional users can use other types)

**skf**: Smoothing kernel function (1 by default, Professional users can use other types)

**sum_den**: (.TRUE. by default, Professional users can use the other type)
- .TRUE.: Use density summation model in the code
- .FALSE.: Use continuity equation

**av_vel**: (.TRUE. by default, Professional users can use the other type)
- .TRUE. : Monaghan treatment on average velocity,
- .FALSE.: No average treatment (C5.99)

**virt_part**: (.FALSE. by default, Professional users can use the other type)
- .TRUE. : Use virtual particle,
- .FALSE.: No virtual particle

**nor_dens**: (.FALSE. by default, Professional users can use the other type)
- .TRUE. : normalized density active
- .FALSE.: No normalized density

**cgra**: Gravity acceleration. (9.8 m/s2 by default)

**dens**: Density of the mixture [Kg/m3]. (It shouldn’t be zero and users should give a value (Range: 1000 to 3000Kg/m3)).

**cmanning**: Voellmy’s coefficient of turbulent viscosity [m/s^2] . (0 by default, User can give a value)

**eros_Coef**: Erosion Coefficient [m-1]. (0 by default, User can give a value)

**nfrict**: Rheological type to calculate basal friction (7 by default, Professional users can use other types)

**Tauy0**: cohesion which is  one of the parameter for Bingham fluids, [N/m^2] (0 by default, User can give a value)

**constK**: It is one of the parameter for Bingham fluids [Pa.s] . (0 by default, User can give a value)

**visco**: viscosity which is  one of the parameter for Bingham fluids [Pa.s] . (0 by default, User can give a value)

**tanfi8**: tangents of the final friction angles [-]. (It shouldn’t be zero. User should give a value)

**hfrict0**: limit for depth of the landslide, to avoid dividing for very small values when obtaining the bottom friction. It should not be zero. (10^(-2) by default, User can give a value)

**c11**: if C11<0 read others 15 constants. (0 by default)

**tanfi0**: Tangents of the initial friction angles [-]. (0. by default, User can give a value)

**Bfact**: vertical consolidation coefficient, Only use if icpwp=1 [s/m2]. (0 by default, User can give a value)

**hrelpw**: The relative width of the basal saturated layer to the total depth (h_w^rel) was assumed to range between 0.25 and 1 [-].  (1 by default, User can give a value)

**Comp**: comparison value for parameter hfrict0. Comp should not be zero. (10^(-2) by default)

**K0**: To active k0 and define the value of this variable is zero if the vertical and horizontal stresses are equal and equivalent to one considering other relationship between the two. In the latter case, the program asks you to enter data for the k (active and passive) that determine this relationship values coefficients. (0 by default, Professional users can give value)

**icpwp**: pore pressure dissipation, Control parameter for interstitial pressure. (0 by default, 1 for activate and 0 for inactive)

**coarse**: when coarse= 1 to build a new coarse mesh for plotting. (0 by default, Professional users can give value)

**chk_pts**: control for output points at which to check solution. (0 by default, Professional users can give value)

**Gid_Mask_SW**: In order to avoid excessive output data, in the filtering vector it will be selected by typing 1 the variables we want to be in the output.

1. **hs**: soil height  (1 for activate and 0 for inactive) 
2. **disp**: displacement   (1 for activate and 0 for inactive) 
3. **v**: velocity   (1 for activate and 0 for inactive) 
4. **Pwb**: Basal pore water pressure   (1 for activate and 0 for inactive) 
5. **eros**: erosion height  (1 for activate and 0 for inactive) 
6. **Z**: Terrain Elevation (1 for activate and 0 for inactive)  
7. **hrel**: Relative height  (1 for activate and 0 for inactive)  
8. **hw**: Water height  (1 for activate and 0 for inactive) 
9. **eta**: porosity  (1 for activate and 0 for inactive)  
10. **hs+hw**: Total height  (1 for activate and 0 for inactive) 
11. **hsat**: Saturated height (1 for activate and 0 for inactive) 
12. **Pw**: Pore-water pressure (1 for activate and 0 for inactive)  

**T_change_to_W**: (1.e+12 by default, Professional users can give value)

# .MASTER.DAT input file:

It is the file giving general inputs to the program. It gives information about the total time the analysis has to run, the analysis time increment and the maximum number of steps after which the analysis stops.

**nline**: indicate the lines of text describing the problem. (1 by default)

**sph**: 1 activates the SPH program (always 1)

**gfl**: 1 activates the FE program (always 0)

**Monte-Carlo**: Type of probabilistic analysis to be performed (0 for the deterministic model, 5 for Monte Carlo and 6 for FOSM)

**problem_type**: General type of problem (always 1): 

**Integ_Alg **: Analysis type (4 for propagtin modeling and 200 for triggering modeling) 

**file.dat**: Problem name (for the .dat file)

**dt** = Time step analysis for reference in output (0.1 by default, Professional users can give value)

**time_end** = Total analysis time for end of computation (1000 by default, Professional users can give value) 

**maxtimestep**: maximum number of timesteps in computation. (10^9 by default, Professional users can give value) 

**print_step**: Number of steps for print (5 by default, Professional users can give value)

**save_step**: Number of steps for save (5 by default, Professional users can give value) 

**plot_step**: Number of steps for plot (5 by default, Professional users can give value) 

**dt_sph**: Time increment which will be used by the code (0.1 by default, Professional users can give value) 

**ic_adapt**: Adaptive control for SPH (always 1)

**time_curves**: If no time curves for histograms will be used, we set ntcurves to zero. (0 by default, Professional users can give value)

**max_pts**: Even if not being used, it is mandatory to provide an integer value. (always 6)

**cases_win**: We can define targets in this line. It is necessary to define the coordinates of one of the angles, along with its length and angle of rotation, for the targets. (0 by default, Professional users can define targets) 

**ic_eros**: Types of bed entrainment laws. (0 by default, Professional users can use different types) 

# .TOP input file:

It is the file that provides topographic data of the problem. 

**ictop**: The type of topography (10 by default, Professional users can other types)

**np**: Number of points (Users should give value)  

**deltx**: The space between the points (Users should give value)

**Topo_x**:  Coordinate of the X-axis

**Topo_y**:  Coordinate of the Y-axis

**Topo_z**: Coordinate of the z-axis

**topo_props**: properties associated to terrain such as basal friction and maximum erodible depth. (0 by default, Professional users can other types)


# .PTS input file:

It is the file that provides the topographic data of the sph nodes with respect to the actual topography of the problem.

**npoin**: indicates the number of points generating the mesh (Users should give value)

**delta**: the preferred distance between them in x and y (Users should give value)

**factshml**: coordinates and the area of influence of every sph node (2 by default, Professional users can give a value)

**Topo_x**:  Coordinate of the X-axis

**Topo_y**:  Coordinate of the Y-axis

**Topo_z**: Coordinate of the z-axis


