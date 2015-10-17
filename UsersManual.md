


---

# Getting started #
## Pre-requisites ##
  * A recent version of Matlab must be installed in the machine meant to run the SHEL software. Currently, matlab version 2009b is being used as the development main machine for the SHEL software, but possibly older versions of matlab will work as well.
  * A Pentium 4 with 1 GB of RAM, or more, and at least 100 MB of available disk space.

## Software access ##
  * The SHEL software package (composed of the source code and some matlab data files) is available in the [downloads section](http://code.google.com/p/shel/downloads/list) or get the latest version from the [subversion repository](http://code.google.com/p/shel/source/checkout).
  * To become an active developer and submit changes to the SHEL project please contact the project owner and a developer account will be provided.
  * Unpack the SHEL software package in a convenient work folder. A few tens of megabytes available should be enough.

## Documentation ##
  * A [technical documentation](TechnicalGuide.md) comes along with the software package in a separate download. It is recommended to download it as well.
  * This wiki.



---

# Quick-start #
  * Download and unzip the latest matlab package available in the [downloads section](http://code.google.com/p/shel/downloads/list) or get the latest version from the [subversion repository](http://code.google.com/p/shel/source/checkout).
  * Open the Matlab environment and set the working folder where the m-files were unpacked.
  * Type 'launch' and and press enter at the matlab prompt.
  * Select a mat-file in the 'initial-conditions' folder and open it.
  * Change the parameters at will.
  * Select two or four-panel view and click the 'show' button.
  * Click the 'run' button to run the model (eventually enable the 'film' option').
  * To export graphics: select a format (eps, png or avi) and a view (all, level, velocity, left or right) from the dropdown menu, then give it a name in the textbox and click  the 'Print' button. If it is the current view only, make sure that the 'film' checkbox is unchecked.



---

# Installation #
  * Simply unpack the SHEL software package in the work-folder of choice.
  * Open the Matlab environment.
  * Set the path to the work folder where the SHEL software was installed.
  * Type _launch_ from the command line and press _enter_. The control panel will show up.



---

# GUI description #

## The control panel ##
![http://content.screencast.com/users/GRiflet/folders/Jing/media/52e8accd-25ef-4e0b-b017-d50f2afcd681/ControlPanelGUIDescrition.png](http://content.screencast.com/users/GRiflet/folders/Jing/media/52e8accd-25ef-4e0b-b017-d50f2afcd681/ControlPanelGUIDescrition.png)

### Physical processes ###
  * **use Coriolis**: activates the Coriolis acceleration term in the momentum equations ie makes the reference frame a rotating reference with the same angular velocity than what is felt at Earths's latitude, given in the _latitude_ parameter, in the _Physics parameters_ panel.
  * **use wind**: activates the wind stress term in the momentum equations. The wind stress is calculated from the wind speed field, given in the _uwind_ and _vwind_ parameters in the _Physics parameters_ panel.
  * **use pressure force**: activates the pressure force due to the elevation gradient. Its mathematical description is provided in the _technical guide_, available in the [download section](http://code.google.com/p/shel/downloads/list).
  * **use bottom drag**: activates the bottom stress, parameterized by the _bottom rugosity_, given in the _physics parameters_ panel.
  * **no-slip**: enforce a no-slip condition, in which the tangential velocity at the lateral walls of the domain is set to zero. This condition re-designs in fact the fluxes masks.

### Physics parameters ###
  * **rho air**: Defines the air density. Enabled if the _use wind_ condition is checked.
  * **rho water**: Defines the water density. Enabled if the _use wind_ condition is checked.
  * **viscosity**: Defines the water turbulent horizontal viscosity. If set to zero, then the  mathematical model solved is equivalent to the Euler equations solution.
  * **gravity**: Earth local vertical gravitational acceleration.
  * **latitude (f)** : Earth latitude of the model domain. It has a correspondence with a rotation rate derived from the Earth's rotation rate. See the _technical guide_ in the [downloads section](http://code.google.com/p/shel/downloads/list) for more details.
  * **bottom rugosity** : Rugosity of the bottom. Will yield the bottom stress value. More details in the _technical guide_ in the [downloads section](http://code.google.com/p/shel/downloads/list).
  * **uwind**, **vwind** : u and v components (along the x-axis and the y-axis respectively) of the wind velocity. Are required to compute the wind stress.
  * **Karman constant**: Adimensional constant required to compute the bottom stress. See the _technical guide_ in the [downloads section](http://code.google.com/p/shel/downloads/list).

### Grid parameters ###
  * **time duration** : Duration of the simulation.
  * **output dt** : Time interval between each graphical output, in the _views panel_.
  * **dt**: Time-step of the momentum numerical scheme.
  * **dx**, **dy**: Space-steps of the regular orthogonal grid, along the x-axis and the y-axis respectively.

### Boundary conditions ###
  * **u-closed**, **v-closed**: when either of these options is checked, the respective sidewalls are considered impermeable, meaning null-flux, and free water elevation. When both these options are checked, the system is similar to a closed domain (closed bathtub) with no open-boundaries.
  * **level**: selects the type of open-boundary condition for the water elevation (ie when the sidewalls are not closed). Two types are available, _Dirichelet_ (or clamped boundary condition) and _Sommerfeld_ (or radiative boundary condition). See the _technical guide_ in the [downloads section](http://code.google.com/p/shel/downloads/list) for more details.
  * **normal**: selects the type of open-boundary condition for the normal velocity to the sidewall (ie when the sidewalls are not closed). Three types are available, _Dirichelet_ (or clamped boundary condition), _Neumann_ (or zero-gradient boundary condition) and _Flather_ (or radiative boundary condition). See the _technical guide_ in the [downloads section](http://code.google.com/p/shel/downloads/list) for more details.
  * **tangential**: selects the type of open-boundary condition for the normal velocity to the sidewall (ie when the sidewalls are not closed). Two types are available, _Dirichelet_ (or clamped boundary condition) and _Sommerfeld_ (or radiative boundary condition). See the _technical guide_ in the [downloads section](http://code.google.com/p/shel/downloads/list) for more details.

### Initial conditions ###
  * **Bump** : generates a gaussian bump in the water elevation, as an initial condition. See the _technical guide_ in the [downloads section](http://code.google.com/p/shel/downloads/list) for the formula of the gaussian elevation.
    * **Lx/2 - x0**, **Ly/2 - y0**: are the distances, in cartesian coordinates (along the x-axis and y-axis), in meters, from the center of the domain to the center of the gaussian bump.
    * **Sx**, **Sy**: are the width of the gaussian bump basis, along the x-axis and along the y-axis.
    * **h**: is the maximum height (in meters) of the gaussian bump elevation.
    * **energy** : is the gaussian bump initial perturbation potential energy. See the _technical guide_ in the [downloads section](http://code.google.com/p/shel/downloads/list) for more details.
    * **Uo** : is the estimated flow velocity in the wake of the circular wave generated by the gaussian bump. See the _technical guide_ in the [downloads section](http://code.google.com/p/shel/downloads/list) for more details.
  * **Taylor col** : is a special initial condition. It computes the initial condition for the geostrophical equilibrium of a pipe flow, flowing over a submarine mount. The user is encouraged to load the _initial-condition_: _taylor-geoIC.mat_ or _taylor-null-IC.mat_, instead of checking this option.
  * **Geostrophic**: _not activated_ in this release. Ignore it.

### Bathymetry ###
  * **island**: Creates a circular-shaped island in the land-mask (where null-fluxes are set).
    * **Lx/2 - x0**, **Ly/2 - y0**: Sets the position of the island's center in distances from the center of the domain to the center of the island, in x-axis and y-axis components.
    * **radius**: Sets the radius of the circular island.
  * **Bathymetry file**: Option that will load the **d** (bathymetry depth) variable from a mat-file, as well as the **eta** (water level) variable (and use it as initial condition).
  * **depth**: Default depth of the domain. Unless some other option is activated, this will be the constant depth of the domain.
  * **M**, **N**: number of grid cells along the x-axis and the y-axis, respectively.
  * **Step-bottom**: Arbitrary splitting of the domain in two regions with constant depth. One region depth is equal to the _depth_ parameter, and the other region's depth is equal to the _step-depth_ parameter.
    * **step-depth**: Depth of the region of the domain with an alternate depth, thus creating a stepping-bottom.

### Tracer ###
  * **tracer**: Activates the advection-diffusion of a passive tracer by the flow field, using an upwind scheme. (To be documented in a future release of the _technical guide_).
    * **K**: Turbulent horizontal diffusion coefficient of the tracer.
    * **Xo**, **Yo**: Initial position of the tracer center of mass relative to the domain center, in cartesian coordinates.
    * **Radius**: Radius of the tracer stain.
    * **Detection treshold**: percentage of the original tracer concentration, above which the grid cell is included when integrating the tracer properties, such as kinetic energy, potential energy, curl, etc...
    * **0 <A,B,C <1**: Displays the _B_ coefficient value (which should be bounded between 0 and 1). The A, B and C coefficients are found in the first and second-order advection-diffusion schemes.
    * **Péclet<1**: Displays the numerical Péclet number of the tracer. See the _technical guide_ in the [downloads section](http://code.google.com/p/shel/downloads/list) for more details.
    * **Umax**: Displays the maximum flow velocity above which the upwind scheme would instabilize. See the _technical guide_ in the [downloads section](http://code.google.com/p/shel/downloads/list) for more details.

### Panels ###
  * **2**: Select the _two-views_ panel. The water elevation view and the flow field view.
  * **4**: Selects the _four-views_ panel. Besides the views from the two-views panel, it adds a _left-property_ view and a _right-property_ view.
  * **Show**: Opens the _views panel_.

<a href='Hidden comment: 
=== Numerical stability ===
=== Caracteristic numbers ===
'></a>

## The views panel ##
![http://content.screencast.com/users/GRiflet/folders/Jing/media/9140125f-7f6d-4d25-aece-639f6e308678/ViewsPanelGUIDescription.png](http://content.screencast.com/users/GRiflet/folders/Jing/media/9140125f-7f6d-4d25-aece-639f6e308678/ViewsPanelGUIDescription.png)

### Time ###
Displays the current time of the simulation in seconds.

### Outputs ###
  * **reset button**: resets the time to 0 and restores the initial conditions of the current configuration.
  * **film checkbox**: check to refresh the views during a running simulation or during a print to avi or to png or eps.
  * **run button**: runs a simulation for the duration set in the _time duration_ textbox, in the _control panel_.
  * **textbox**: must be edited to define the name of the sub-folder, containing the images and video results, under the _image-results_ folder.
  * **dropdown menu**: selects the output format (**avi**, **png** and **eps**) and the views to export (all the _views panel_, the _level view_, the _velocity view_, the _left property view_ and the _right property view_). The _film checkbox_ wil automatically be checked when selecting _avi_ format, and shouldn't be changed while the _avi_ format is chosen.
  * **print button**: Will export the current view(s) in the chosen format (selected via the _dropdown menu_) and store them within the _images-results_ folder (in the sub-folder defined in the _textbox_). If the _avi_ format is chosen (or if the _film checkbox_ is activated) then it will run a simulation (like the _run button_ does) and it will export every frame output in the view(s).

### View sub-panel ###
  * **z-scale**: slider that sets _n_, the power of ten of the maximum and minimum value of the z-axis scale, for the 3D views. (-10<sup>n</sup> 10<sup>n</sup>) for each sub-panel. Present in all sub-panels, except in the _velocity_ sub-panel.
  * **dropdown menu**: Dropdown menu that selects a property from the _property list_ to view. Present in the _left property_ and _right property_ sub-panels only.

### Views ###
  * **level**: Top-left sub-panel. Displays the water level property (_eta_) only, in 3D view.
  * **velocity**: Top-right sub-panel. Displays the vector field of the water flow _(Hu, Hv)_, in a 2D vector plot.
  * **Left property**, **Right property**: Bottom-left and bottom-right sub-panels, respectively. Can display any property from the _property list_.

### Properties view-types ###
Each property of the  _property list_ will assume one of the following types of views:
  * **vector plot**: Vector plot of a 2D vector field.
  * **3D surface**: 3D surface, rendered with a colorbar, by default.
  * **color map**: 2D color map, rendered with a colorbar, by default.
  * **time-plot**: Time-serie plot.
  * **cross-section plot**: plot of a cross-section of the domain.
  * **contour map**: 2D contour map, rendered with a legend of the contour lines, by default.

### Property list ###
This is the exhaustive list of the available properties to render in the sub-panels of the _views panel_.
  1. **bathymetry** : is a _color map_ that shows the bathymetric depth (stored in the _`d`_ variable).
  1. **eta** : is a _3D surface_ that shows the water elevation (stored in the _`eta`_ variable).
  1. **(u, v)** : is a _vector plot_ that shows the water velocity (stored in the _`u_t`_ and _`v_t`_ variables).
  1. **(hu, hv)** : is a _vector plot_ that shows the water flow (based on the product of the water depth, _`H`_, times the _`u_t`_ and _`v_t`_ variables).
  1. **land mask** : is a _color map_ that shows the land mask (stored in the _`mask`_ variable).
  1. **U mask** : is a _color map_ that shows the land mask (stored in the _`mask_u`_ variable).
  1. **V mask** : is a _color map_ that shows the land mask (stored in the _`mask_v`_ variable).
  1. **W mask** : is a _color map_ that shows the land mask (stored in the _`mask_z`_ variable).
  1. **u** : is a _3D surface_ that shows the water elevation (stored in the _`u_T`_ variable).
  1. **v** : is a _3D surface_ that shows the water elevation (stored in the _`v_t`_ variable).
  1. **velocity modulus** : is a _3D surface_ that shows the water elevation (computed with the _`u`_ and _`v`_ variables).
  1. **ke** : is a _3D surface_ that shows the water elevation (stored in the _`ke`_ variable).
  1. **pe** : is a _3D surface_ that shows the water elevation (stored in the _`pe`_ variable).
  1. **e** : is a _3D surface_ that shows the water elevation (sum of _`ke`_ and _`pe`_ variables).
  1. **eke** : is a _3D surface_ that shows the water elevation (stored in the _`eke`_ variable).
  1. **curl** : is a _3D surface_ that shows the water elevation (stored in the _`curl_t`_ variable).
  1. **enstrophy** : is a _3D surface_ that shows the water elevation (stored in the _`enstrophy_t`_ variable).
  1. **squared shear** : is a _3D surface_ that shows the water elevation (stored in the _`sqshearrate_t`_ variable).
  1. **squared strech** : is a _3D surface_ that shows the water elevation (stored in the _`sqstrechrate_t`_ variable).
  1. **squared strain** : is a _3D surface_ that shows the water elevation (stored in the _`sqstrain_t`_ variable).
  1. **Okubo-Weiss** : is a _3D surface_ that shows the water elevation (stored in the _`okuboweiss_t`_ variable).
  1. **global volume** : is a _time plot_ that shows the perturbation volume of the domain (stored in the variable _`volume`_).
  1. **global momentum** : is a _time plot_ that shows the domain integrated u and v component of the velocity field (stored in the variables _`iMomentumU`_ and _`iMomentumV`_).
  1. **global e-ke-pe-eke** : is a _time plot_ that shows the domain integrated energy (stored in the variables _`iKe`_, _`iPe`_ and _`iEke`_).
  1. **global vorticity** : is a _time plot_ that shows domain integrated relative vorticity (stored in the variable _`iVort`_).
  1. **global enstrophy** : is a _time plot_ that shows domain integrated relative enstrophy (stored in the variable _`iEnst`_).
  1. **global squared shear** : is a _time plot_ that shows the domain integrated squared shear (stored in the variable _`iSqShear`_).
  1. **global squared strech** :  is a _time plot_ that shows the domain integrated squared strech (stored in the variable _`iSqStrech`_).
  1. **global squared strain** : is a _time plot_ that shows the domain integrated squared strain (stored in the variable _`iSqStrain`_).
  1. **global Okubo-Weiss** : is a _time plot_ that shows the domain integrated squared strain, enstrophy and Okubo-Weiss scalar (stored in the variables _`iSqStrain`_, _`-iEnst`_ and _`iOWeiss`_).
  1. **Section X-Axis Water level**: is a _cross-section plot_ along the x-axis that shows the water elevation profile contour (stored in the variable _`eta`_).
  1. **Okubo-Weiss (contour)**: is a _contour map_ of the Okubo-Weiss scalar, where positive values are represented by dashed contour lines, and negative values are represented by a solid contour line. Zero values are represented by a thicker solid contour line.
  1. **Tracer** :  is a _color map_ that shows the tracer concentration (stored in the _`Tr`_ variable).



---

# How to #
## Load a configuration ##
  * From the control panel click the _load_ button, or from the matlab command line prompt, type _launch_ and press _enter_. It will load the _control panel_.
  * Then choose a saved configuration file (with a .mat extension). By default, a few configuration samples come built-in with the SHEL software package in the _initial conditions_ folder. Click on _open_.
![http://content.screencast.com/users/GRiflet/folders/Jing/media/c5628f48-5c1a-478b-bbce-b561983040a2/controlpanel_openfromload.png](http://content.screencast.com/users/GRiflet/folders/Jing/media/c5628f48-5c1a-478b-bbce-b561983040a2/controlpanel_openfromload.png)

## Visualize a configuration ##
  * To visualize a loaded configuration, select one of the _two_ or the _four_ radio-button (enabling the two-panel or the four-panel view), then click on button _show_.
![http://content.screencast.com/users/GRiflet/folders/Jing/media/be7de457-a970-4a61-b3de-4972e898a7bc/controlpanel-show.png](http://content.screencast.com/users/GRiflet/folders/Jing/media/be7de457-a970-4a61-b3de-4972e898a7bc/controlpanel-show.png)
  * In the _view panels_, the user can run the model, select other views (if on the four-view panel), select the scale, reset the model to the initial instant and export the results in image files or avi movie file.

## Run the model ##
  * Click on the _reset_ button at any time to reload the initial instant of the loaded configuration.
  * Click on the _run_ button to compute all the time instants (intervalled by the _dt_ parameter set in the _dt_ textbox in the _control panel_). If the _film_ checkbox is checked, then every _output dt_ (a parameter set in the _output dt_ textbox, in the _control panel_), the simulation results will be printed onscreen in the two-panel or the four-panel views. Graphical output may slow down considerably the simulation and depends on the machine hardware (which gpu, which cpu, ...). Thus, to speed up the simulation duration, simply uncheck the _film_ checkbox at any time.
  * The simulation will last a simulated duration equal to the value set in the _time duration_ textbox, in the _control panel_. Then, it will graphically display the last computed instant in the _view panels_, whether the _film_ checkbox is or is not checked.
![http://content.screencast.com/users/GRiflet/folders/Jing/media/3d6c1b34-bdb5-423f-be0c-90391b2dba8d/viewpanels-run.png](http://content.screencast.com/users/GRiflet/folders/Jing/media/3d6c1b34-bdb5-423f-be0c-90391b2dba8d/viewpanels-run.png)

## Export graphics and movies ##
To help users to produce paper-standard captions, it is possible to export results in both vectorial and raster formats, respectively **eps** and **png**. It is also possible to export in **avi** movie format.
  1. Insert a name near the _print_ button textbox. It will create a new folder in the _image-results_ work folder, with that name, where it will place all subsequent images and movies.
  1. Select which variables to display from the drop down menus, _left property_ and _right property_ (in the four-panel view only).
  1. Set the scale of each view with the respective _scale slider_.
  1. Rotate and adjust the zooming of the camera, for each view, by hitting the _rotate 3D_ icon in the top _menubar_.
  1. Select which view and which format to print from the dropdown menu, near the _reset_ button. Each view of the two(_level_, _velocity_), or four (_left_, _right_), available views is selectable. Alternatively, the whole views are selectable (_all_). The available output formats are **png** (raster), **eps** (vectorial) or **avi** (movie).
![http://content.screencast.com/users/GRiflet/folders/Jing/media/3e661745-ed14-4ac1-ad56-0e0fa4a1bf25/viewpanels-export1.png](http://content.screencast.com/users/GRiflet/folders/Jing/media/3e661745-ed14-4ac1-ad56-0e0fa4a1bf25/viewpanels-export1.png)
  1. Optionally, if the _film_ checkbox is checked then a simulation will be continued and will last for the whole _time duration_, set in the _control panel_, and each _output dt_ will be exported in the chosen media format. Note that for **avi**, the _film_ checkbox will automatically be set and should not be unchecked during the frame scraping. If the _film_ checkbox is unchecked, then only what is currently displayed on the view panels will be exported.
  1. To export the graphics click on the _print_ button then wait while the process takes place. It may take some time depending on the format (**png** are the slowest), the type of view and if the _film_ checkbox is set.
![http://content.screencast.com/users/GRiflet/folders/Jing/media/14f3b79d-d6a9-4aab-87fd-ce1602275992/viewpanels-export2.png](http://content.screencast.com/users/GRiflet/folders/Jing/media/14f3b79d-d6a9-4aab-87fd-ce1602275992/viewpanels-export2.png)
  1. After the process of exporting the graphics took place, the newly saved graphics are found in the _image-results_ folder, in the sub-folder with the same name given in the textbox near the _print_ button.
  1. Remember:
    1. because it is very important for the user to reproduce the same results (and the same graphics) at a later stage, both the _control panel_ configuration is saved, along with the _view panels_ settings, in a mat-file with the same name and in the same sub-folder. This allows the user to re-open the exact same configuration, and view settings, during another matlab session.
    1. Additionally, the one-dimensional integrated diagnostics time-series (integrated energy, vorticity, etc...) are also kept in a mat-file with the same name, added with the suffix _global_. This allows the user to re-plot the time-series, later, accordingly to the editorial standards of a given journal.
![http://content.screencast.com/users/GRiflet/folders/Jing/media/584480c8-07b8-4cb5-8af6-0afafa326b69/viewpanels-resultsfolder.png](http://content.screencast.com/users/GRiflet/folders/Jing/media/584480c8-07b8-4cb5-8af6-0afafa326b69/viewpanels-resultsfolder.png)

## Edit and save a new configuration ##
To save a new configuration (from the _control panel_):
  1. Edit the parameters in the _control panel_ GUI.
  1. Eventually (advanced users only), use the matlab debug tool to manually edit the bathymetry matrix and the initial conditions matrices.
  1. Click the _save_ button in the _control panel_.
  1. From the new dialog box, select the folder _initial-conditions_.
  1. In the textbox, insert a name for the new configuration.
  1. Click the _save_ button of the dialog box.
![http://content.screencast.com/users/GRiflet/folders/Jing/media/33a92fb2-6a17-4597-9b37-ec46d89629c6/controlpanel-saveconf.png](http://content.screencast.com/users/GRiflet/folders/Jing/media/33a92fb2-6a17-4597-9b37-ec46d89629c6/controlpanel-saveconf.png)

Alternatively, it was seen, in the section _[how to export graphics and movies](#Export_graphics_and_movies.md)_, how to save the configuration, in the _control panel_, along with the settings, in the _views panel_, all in the same folder than of the graphics.

<a href='Hidden comment: 
= Advanced: How to create a generic ...=
To be done...
== matricial variable ==
== bathymetry ==
== land-mask ==
== initial-condition ==
== tracer ==
'></a>



---

# To learn more about SHEL #
Either explore the GUI, either contact the author: `guillaume.riflet at gmail.com`.