

---

# Design philosophy #
Here are some development guidelines:
  * **Devise coherent global variables that have a specific role in numerical modeling of fluid mechanics.** For example, `eta` is the water elevation, `d` is bathymetry depth and `H` is the water depth (sum of `eta` and `d`). These global variables can be invoked in any part of the code and, thus, their list must be constantly updated and synchronized between this wiki and the code.
  * **Separate functions, or replaceable parts of the numerical engine, in different files, but keep the number of files as low as possible.** For example, the temporal numerical scheme is contained in the file `ComputeLeapfrog.m`, the spatial numerical scheme is contained in the file `ComputeSpaceU_CS.m`, and the open boundaries schemes are contained in the file `computeOB_U.m`. That way, it is easy for other developers to replace the current functional file with another. It makes it also easy to upgrade to a system where the user can choose between different types of schemes for the same function; like for example, to change the `ComputeLeapfrog.m` file with a `ComputeMcCormack.m`.
  * **Use the Matlab index notation to perform matricial computations. `do loops` are forbidded!** The `profiler on` tool showed a 10x improvement in computation speed when `do loops` are religiously avoided.
  * **Make the graphical interface as simple, friendly, flexible and bullet-proof as possible.** The goal is to minimize the number of clicks to achieve an action, but also to opt-out possible pathways that would lead the user to a dead-end. For example, once the _tracer_ checkbox is deactivated in the _control panel_, then all the options related becomed disabled (greyed out). This is powerful visual feedback to the user about his choices left.
  * **Infinite configurations are possible and yet, a past successful configuration must be always reloadable.** This means that the state of the global variables that define a configuration must be easily saved and later reloaded. This includes matricial fields like the bathymetry, `d`, or the velocity fields, `u` and `v`.
  * **Whenever an image is exported, its configuration file must be saved along with it.** This is to allow the user to later reproduce the image and eventually change the legend, or some of its parameters. It is better than saving the whole fig file for each image.

## Futur design improvements ##
This is mostly a wish-list of futur improvements in the design and structure of the SHEL interface:
  * **Integrate the _control panel_ in the _views panel_ as a _menubar_ with child menus**. Users would get bigger a sense of familiarity with the SHEL interface.
  * **Load a preview snapshot of saved configurations**, added of a short description by the creator.
  * **Find a way to** centralize, **distribute** and promote **SHEL configuration files**.
  * **Allow users to choose from different numerical schemes from dropdown menus.** Make it easy from the developers side to integrate new numerical schemes in the SHEL.
  * **See if the code is faster by avoiding to use the `global` directive**, while keeping in mind the same global variables philosophy. Perhaps by pushing structures in the function calls...

---

# The Arakawa C staggered grids matricial design #
The SHEL is about 2D layers of water properties, be it elevation, u-component of velocity, v-component of velocity, energy or enstrophy. Thus, most critical variables will be stored in matrices of specific sizes, representing a specific mesh, which is regular, orthogonal, cartesian and ... staggered. This mesh is well known in numerical modeling literature as the staggered Arakawa C grid. There are four types of orthogonal, rectangular, regularly spaced grids: the **T-cells**, the **U-cells**, the **V-cells** and the **W-cells**. Each are natural discretized representations of the computed variables. For example the u-component of velocity are best represented in the U-cells, however, the water elevation is best represented in a T-cell. They are sketched as follows:

```
'%   --j y
%   |           --- V ---           V ------- V         U -- T -- U
%   i          |         |          |         |         |         |
%   x          U    T    U          T    U    T         |    V    |
%              |         |          |         |         |         |
%               --- V ---           V ------- V         U -- T -- U
%                 T-cell               U-cell              V-cell
%                  M x N               M+1 x N             M x N+1
%       T(eta, H, tr, d, f, mask,   U(u, mask_u         V(v, mask_v
%           x, y, wind and bottom      x_u, y_u)           x_v, y_v)
%            stress, curl(v,u,0))
%
%    --j y
%   |           --- U ---  
%   i          |         | 
%   x          V    W    V 
%              |         | 
%               --- U ---  
%                 T-cell
%                M+1 x N+1
%       W( curl(u,v,0) )
%        (deformation rate = curl(-u,v,0)'
```

  * **T-cells**:
    * size: MxN.
    * Variables: water elevation, bathymetry depth, water column depth, land mask.
  * **U-cells**:
    * size (M+1)xN.
    * variables: u-component of velocity, flux mask along u-axis.
  * **V-cells**:
    * size Mx(N+1).
    * variables: v-component of velocity, flux mask along v-axis.
  * **W-cells**:
    * size (M+1)x(N+1).
    * variables: relative vorticity (curl).

Bear in mind that it is often required to convert a variable from one grid type to another. However, converting one variable from one grid type to another (say, from the W-cell to the T-cell), requires an interpolation algorithm. In the SHEL code, the simplest linear interpolation algorithm are used. They perform an average with two or four grid points, depending on the source and target grids of the interpolation.

**Below is the interpolation algorithm from the W-cells grid to the T-cells grid:
```
'%%File: ComputeDiagnostics.m'
function tgrid = interpolFtoT(fgrid, mask, M, N)
tgrid = .25 * mask .* ( ...
              fgrid(1:M,1:N) ...
            + fgrid(2:M+1,1:N) ...
            + fgrid(2:M+1,2:N+1) ...
            + fgrid(1:M,2:N+1) ...
            );
```**

**Below is the interpolation algorithm from the V-cells grid to the U-cells grid and from the T-cells grid to the U-cells grid:
```
'%%File ComputeSpaceU_CS.m'
function av = fouraverage_u(v, M, N) 
    av =  .25 * ( ... 
            v(2:M,1:N) + ... 
            v(2:M,2:N+1) + ... 
            v(1:M-1,1:N) + ... 
            v(1:M-1,2:N+1) ... 
           ); 
         
function av = twoaverage_u( H, M, N)
    av = .5 * ( ...
            H(2:M,2:N-1) + H(1:M-1,2:N-1) ...
         );
        
```**

Please refer to the _technical guide_ in the [downloads section](http://code.google.com/p/shel/downloads/list) for more details.


---

# Resources listing #
Here's an exhaustive description of the SHEL inner resources, such as folders, files, functions and global variables.

## Folders ##
  * `root` folder: contains the SHEL matlab scripts.
  * `initial-conditions`: contains SHEL configurations in the matlab data file format (.mat). See the TemplatesDescription for more details
  * `images-results`: if not present, the folder is created after the user exports images from the SHEL views panel. It contains a set of sub-folders, each containing images, avi movies generated by the user, along with the SHEL configuration file and the timeseries data files.
  * `bathymetries`: contains files containing bathymetry (`d`), level (`eta`) and velocity (`u`, `v`) variables stored in the matlab data file format (.mat).

## Files ##
  * [launch.m](http://code.google.com/p/shel/source/browse/trunk/launch.m): fires up the _control panel_ and its _resetbutton_ callback.
  * [ControlPanel2D](http://code.google.com/p/shel/source/browse/trunk/ControlPanel2D.m): codes the full _control panel_ interface.
    * [About.m](http://code.google.com/p/shel/source/browse/trunk/about.m): Displays the GPL license agreement.
    * [visualize\_v2.m](http://code.google.com/p/shel/source/browse/trunk/visualize_v2.m): codes the full _two-views panel_ interface.
    * [visualize\_v4.m](http://code.google.com/p/shel/source/browse/trunk/visualize_v4.m): codes the full _four-views panel_ interface.
      * [initialconditions.m](http://code.google.com/p/shel/source/browse/trunk/initialconditions.m): codes the initial conditions for the global matricial variables.
      * [ComputeModel.m](http://code.google.com/p/shel/source/browse/trunk/ComputeModel.m): codes the main loop in time of a SHEL simulation.
        1. [ComputeLeapfrog.m](http://code.google.com/p/shel/source/browse/trunk/ComputeLeapfrog.m): codes the momentum and water level after a single iteration in time.
          * [ComputeContinuity.m](http://code.google.com/p/shel/source/browse/trunk/ComputeContinuity.m): codes the water level increment during a single time-iteration.
          * [ComputeSpaceU\_CS.m](http://code.google.com/p/shel/source/browse/trunk/ComputeSpaceU_CS.m): codes the momentum increments during a single time-iteration-
          * [ComputeTimeU\_FT.m](http://code.google.com/p/shel/source/browse/trunk/ComputeTimeU_FT.m): codes the new momentum after the single time-iteration for the inner-domain.
          * [computeOB\_U.m](http://code.google.com/p/shel/source/browse/trunk/computeOB_U.m): codes open-boundary conditions.
        1. [ComputeDiagnostics.m](http://code.google.com/p/shel/source/browse/trunk/ComputeDiagnostics.m): codes the computation of the diagnostic quantities.
        1. [ComputeTracer\_FT.m](http://code.google.com/p/shel/source/browse/trunk/ComputeTracer_FT.m): codes the tracer after a single time-iteration.
        1. [plotmodel.m](http://code.google.com/p/shel/source/browse/trunk/plotmodel.m): codes the views graphical rendering.
        1. [printit.m](http://code.google.com/p/shel/source/browse/trunk/printit.m): codes the exporting to png, eps or avi formats.

## Functions ##
This list of function is yet unexhaustive. Only the main functions of the SHEL numerical modeling core are described. The interface functions from the Matlab GUI are _not_ described.
  * `launch`: fires up the _control panel_ and its _resetbutton_ callback.
  * `controlpanel2D`: GUI interface of the _control panel_.
    * `About`:
    * `visualize_v2`: GUI interface of the _two-views panel_.
    * `visualize_v4`: GUI interface of the _four-views panel_.
      * `initialconditions`:
        1. `updatetime`:
        1. `makecoordinates`:
        1. `fillfields`:
      * `ComputeModel`:
        1. `ComputeLeapfrog`:
          * `RHSeta = ComputeContinuity(H, eta, u, v, mask)`:
            1. `RHS = ComputeContinuityU(H,u,mask,dx)`:
          * `RHSu = ComputeSpaceU_CS(H, H_old, eta, u, u_old, v, v_old, mask_u, dx, dy, signcoriolis)`:
            1. `av = fouraverage_u(v, M, N) `:
            1. `av = twoaverage_u( H, M, N) `:
            1. `tb_L = bottomstress( u_L, v_L, H_L, lb, karman) `:
            1. `u_L = windstress(u_L, v_L, rho0, rho_air) `:
          * `u_new = ComputeTimeU_FT(u_old, u_new, RHSu, H_old, H_new, mask_u, dt)`:
          * `[eta_new, H_new, u_new, v_new] = computeOB_U( eta_new, eta_old, H_new, H_old, d, u_new, mask_u, v_new, v_old, mask_v,dt, dx)`:
        1. `ComputeDiagnostics(l)`:
          * `computeCurl_w`:
          * `computeShearRate_w`:
          * `computePotentialVorticity_w`:
          * `tgrid = interpolFtoT(fgrid, mask, M, N)`:
        1. `ComputeTracer_FT(Tr, K)`:
          * `RHSt = ComputeSpaceT_UP_U(H, u, mask_u, Tr, K_L, dx_L)`:
          * `RHSt = ComputeFaceFluxT_UP(Hm, TrUp, TrUm, TrDif, dx_L, K_L)`:
        1. `plotmodel(handles)`:
          * `plotproperty(prop,CLim,Hn)`:
            1. `SetXY(Hn)`:
            1. `SetZ(CLim, label, Hn, az, el)`:
        1. `printit(frame_L, handles, choice)`:

## Global variables ##

### loadable and save-able ###
Excerpt taken from the [ControlPanel2D.m](http://code.google.com/p/shel/source/browse/trunk/ControlPanel2D.m) file:
```
'%Physical processes'
global coriolis; '%use coriolis force? Logical.'
global bottom; '%use bottom stress? Logical.'
global wind; '%use wind stress? Logical.'
global pressure; '%use pressure gradient force? Logical.'
global noslip; '%use no-slip? Logical.'

'%Physics parameters'
global rho_air; '%air density. Real.'
global rho0; '%water density. Real.'
global Nu; '%water turbulent horizontal viscosity. Real.'
global g; '%gravitational pull. Real.'
global f; '%Coriolis frequency. Real.'
global lb; '%Bottom rugosity length. Real.'
global uwind; '%x-axis wind speed component. Real.'
global vwind; '%y-axis wind speed component. Real.'
global karman; '%Von Karman constant. Real.'

'%Grid parameters'
global duration; '%time duration. Real'
global outputL; '%number of iterations for output. Hidden. Integer.'
global outputdt; '%time interval between graphical outputs. Real.'
global dt; '%time step. Real.'
global dx; '%x-axis grid step. Real.'
global dy; '%y-axis grid step. Real.'

global gama; '%gama coefficient for Asselin-Robert filter. Hidden. Real.'

'%Bathymetry'
global d0; '%Constant depth. Real.'
global M; '% x-axis grid size. Integer.'
global N; '% y-axis grid size. Integer'
global L; '% Characteristic length. Real.'
 '%Bathymetry file'
global loadbathymetry; '%bathymetry test. Logical.'
global bathymetryfile; '%bathymetry file. String.'
global d; '%The current bathymetry. Matrix.'
global eta0; '% Initial water elevation. Matrix.'
 '%step bottom'
global step; '%Step test. Logical.'
global d0_step; '%Step depth. Real.'
 '%island'
global tc_isla; '%Use island? Logical.'
global isla_x0; '%island x-position. Real.'
global isla_y0; '%island y-position. Real.'
global isla_R; '%island radius. Real.'

'%Initial conditions'
global tc_taylor; '%Use taylor column initialization? Logical.'
global tc_geostrophic; '%user geostrophic initial fields? Useless. Logical.'
 '%bump'
global tc_bump; '%Use gaussian elevation? Logical.'
global bump_x0; '%gaussian center x-position. Real.'
global bump_y0; '%gaussian center y-position. Real.'
global bump_sx; '%gaussian width along x-axis. Real.'
global bump_sy; '%gaussian width along y-axis. Real.'
global bump_d0; '%gaussian max height. Real.'

'%tracer'
global tracer; '%Use tracer? Logical.'
global K; '%tracer horizontal turbulent diffusivity. Real.'
global TrXo; '%tracer center x-position. Real.'
global TrYo; '%tracer center y-position. Real.'
global TrR; '%tracer radius. Real.'
global Treshold; '%tracer min concentration, above which integration is '
                 '%made in diagnostics. Real.'

'%boundary conditions'
global radiatelevel; '%Use level radiation? Logical. '
global flather; '%Use flather radiation? Logical.'
global neumann; '%Use neumann condition? Logical.'
global radiatetan; '%Use radiation for tangential velocity? Logical.'
global u_closed; '%Close the x-axis boundary? Logical.'
global v_closed; '%Close the y-axis boundary? Logical.'

'%Panels outputs & properties'
global npanels; '%Two or four panels? Integer.'
global leftprop;  '%Property number from property list. For the left panel. Integer.'
global rightprop; '%Property number from property list. For the right panel. Integer.'
 '%CLim coef A: [-1eA 1eA]'
global CLimEta; '% m [-10. 2.] Water level Z-scale limits. Real vector.'
global CLimV; '% m/s [0. 100.] velocity Z-scale limits. Real vector.'
global CLimLeft; '% J [-10. 10.] left property z-scale limits. Real vector.'
global CLimRight; '% J [-10. 10.] right property z-scale limits. Real vector.'
 '%print to file'
global film;  '%Do we animate the simulation and render it on screen? Logical.'
global printG; '%Do we export the screen rendering to a file? Logical.'
global myfile; '%filename to build sub-folder and image filenames. String.'
global optprint; '%Which export format and which property? Integer.'
global movie; '%Is it the avi format we want to export? Logical.'
global frame; '%frame counter of the graphical exports. integer.'
```

### Initializable global matrices ###
_Please refer to the [technical guide](TechnicalGuide.md) for the mathematical formulation and the numerical algorithms for the computation of prognostic and diagnostic variables. In particular, the diagnostic variables computation is coded in the [computeDiagnostics.m](http://code.google.com/p/shel/source/browse/trunk/ComputeDiagnostics.m) file._

Excerpts taken from [initialconditions.m](http://code.google.com/p/shel/source/browse/trunk/initialconditions.m).
```
'%function initialconditions'
global time; '% time variable. Real.'
global timetr; '% tracer time varialbe. Real.'

'%function makecoordinates'
global x; '%x-component of the T-grid coordinates. Real matrix.'
global y; '%y-component of the T-grid coordinates. Real matrix.'

global x_u; '%x-component of the U-grid coordinates. Real matrix.'
global y_u; '%y-component of the U-grid coordinates. Real matrix.'

global x_v; '%x-component of the V-grid coordinates. Real matrix.'
global y_v; '%y-component of the V-grid coordinates. Real matrix.'

global x_w; '%x-component of the W-grid coordinates. Real matrix.'
global y_w; '%y-component of the W-grid coordinates. Real matrix.'

'%function fillfields'
'%T-cells'
global Tr; '% tracer concentration.'
global eta0; '% initial water elevation.'
global eta_old; '%water elevation at instant t-dt.'
global eta; '%water elevation at instant t.'
global eta_new; '%water elevation at instant t+dt.'
global H_old;' %water column height at instant t-dt.'
global H;' %water column height at instant t.'
global H_new; '%water column height at instant t+dt.'
global d; '%bathymetry depth. constant in time.'
global mask; '%land mask.'
global masknan;' %not-a-number mask.'
global u_t;' %u component of velocity.'
global v_t; '%v component of velocity.'
 '%Diagnostic quantities'
global ke; '%kinetic energy.'
global pe; '%potential energy.'
global eke; '%turbulent kinetic energy.'
global gradux;' %u gradient along the x-coordinate.'
global graduy; '%u gradient along the y-coordinate.'
global gradvx; '%v gradient along the x-coordinate.'
global gradvy;' %v gradient along the y-coordinate.'
global curl_t;' % (v_x - u_y) The curl of velocity'
global potvorticity_t; '% the potential vorticity'
global strechrate_t;' % strech rate. (u_x - v_y) %Check Arakawa 1966'
global shearrate_t; '% shear rate. (v_x + u_y) %Check Arakawa 1966'
global sqstrechrate_t; '%the square of the strech rate.'
global sqshearrate_t; '%the square of the shear rate.'
global enstrophy_t; '% enstrophy (0.5 * curl * curl)'
global sqstrain_t;' % the square of the strain rate. '
               ' % 0.5 * ( shear * shear + strech * strech )'
global divergence_t; '%the horizontal divergence.'
global sqdivergence_t;' %the square of the horizontal divergence.'
global okuboweiss_t; '% the okubo-weiss scalar (sqstrain - enstrophy) '
                    '%(Check Arakawa1966 and Weiss1981)'

'%U-cells'
global u_a; '%time averaged u component of velocity.'
global u_old; '%u component of velocity at instant t-dt.'
global u; '%u component of velocity at instant t.'
global u_new; '%u component of velocity at instant t+dt.'
global mask_u; '%flux mask of the u-component of velocity.'

'%V-cells'
global v_a;' %time averaged v component of velocity.'
global v_old; '%v component of velocity at instant t-dt.'
global v; '%v component of velocity at instant t.'
global v_new; '%v component of velocity at instant t+dt.'
global mask_v; '%flux mask of the v-component of velocity.'

'%Z-cells'
global curl_w; '%The curl of velocity (v_x - u_y).'
global potvorticity_w; '%Rossbys potential vorticity.'
global shearrate_w;' % shear rate. (v_x + u_y) %Check Arakawa 1966.'
global sqshearrate_w; '% the square of the shear rate.'
global enstrophy_w; '% enstrophy (0.5 * curl * curl).'
global okuboweiss_w; '% the okubo-weiss scalar (sqstrain - enstrophy) '
                   ' %(Check Arakawa1966 and Weiss1981).'

'%time-dependent global properties (domain integrated)'
global volume; '%global volume.'
global iKe; '% kinetic energy.'
global iPe; '% potential energy.'
global vtime; '% time vector.'
global iVort; '% relative vorticity.'
global iEnst; '% enstrophy.'
global iSqStrech; '%square of the strech rate.'
global iSqShear;' %square of the shear rate.'
global iSqStrain; '%square of the strain rate.'
global iOWeiss; '%okubo-weiss parameter.'
global iMomentumU; '%the u-component of momentum. '
global iMomentumV;' %the v-component of momentum.'
global iEke; '% turbulent kinetic energy.'
```


---

# To learn more about SHEL development #
Please, contact the author: `guillaume.riflet at gmail.com`.