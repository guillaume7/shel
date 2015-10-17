

# Motivation #
Often, the user needs to compute a new diagnostic quantity, one that is not available on the current [property list](UsersManual#Property_list.md). Hence, this tutorial serves to guide the user to become a developer, and to help him create new diagnostic quantities in the SHEL. The developer is also encouraged to submit new quantities to this googlecode project.

**In this tutorial, a dummy property, called `dummy1`, will be created.**

# Which grid type? #
**Determine which is the best suited [grid type](DevelopersManual#The_Arakawa_C_staggered_grids_matricial_design.md)**.

All properties are, eventually, interpolated on a [T-cell](DevelopersManual#The_Arakawa_C_staggered_grids_matricial_design.md), for graphical rendering purposes. However this does not mean they are all computed on T-cells. This is the first important design step. Is this quantity naturally computed on T, U, V or W-cell grid?

For example, the `curl` is a quantity naturally computed on the W-cell (see the [technical guide](TechnicalGuide.md) for more details on the which grid cells are computed the properties). Thus, the global variable `curl_w` was [computed on a W-cell and was later interpolated on a T-cell](http://code.google.com/p/shel/source/browse/trunk/ComputeDiagnostics.m#233), to the global variable `curl_t`.

**In this case, `dummy1` is computed on a W-cell.**

# Create the diagnostic variable #

Edit the following files.

## initialconditions.m ##

Create the global variable `dummy1_t`, at [the end of the diagnostic quantities in the T-cells code zone](http://code.google.com/p/shel/source/browse/trunk/initialconditions.m#196),
```
global sqdivergence_t; '%the square of the horizontal divergence.'
global okuboweiss_t; '% the okubo-weiss scalar (sqstrain - enstrophy) '
                    '%(Check Arakawa1966 and Weiss1981)'
global dummy1_t; '% always comment what this quantity is!'
```

and create the global variable `dummy1_w`. [at the end of the W-cells (or Z-cells )code zone](http://code.google.com/p/shel/source/browse/trunk/initialconditions.m#232),
```
global okuboweiss_w; '% the okubo-weiss scalar (sqstrain - enstrophy)' 
                    '%(Check Arakawa1966 and Weiss1981)'
global dummy1_w; '% always comment what this quantity is!'
```

Allocate memory for the new global variable `dummy1_t`, [at the end of the T-cells allocation](http://code.google.com/p/shel/source/browse/trunk/initialconditions.m#442):
```
sqdivergence_t = zeros(M,N);
dummy1_t = zeros(M,N);
```

and allocate memory for the new global variable `dummy1_w`, [at the end of the W-cells allocation (Z-cells)](http://code.google.com/p/shel/source/browse/trunk/initialconditions.m#416):
```
okuboweiss_w = zeros(M+1,N+1);
dummy1_w = zeros(M+1,N+1);
```

## computeDiagnostics.m ##

Create the function that computes the `dummy1_w` global variable, at the end, [right before the `interpolFtoT` function](http://code.google.com/p/shel/source/browse/trunk/ComputeDiagnostics.m#383).

```
'%%%% New code begins here %%%%'
function computeDummy1_w
  global u;
  global v;
  global dummy1_w; '%W-cell global variable'
  global dummy1_t; '%T-cell global variable'
  global mask; '%T-cell land mask'

  '%Algorithm of dummy1'
  dummy1_w(2:M,2:N) = u(1:M-1,2:N) + u(2:M,2:N) + v(2:M,1:N-1) + v(2:M,2:N);

  '%Interpolation of dummy1 to a T-cell'
  dummy1_t = interpolFtoT(dummy_w, mask, M, N)
'%%%% New code ends here %%%%'

function tgrid = interpolFtoT(fgrid, mask, M, N)
'%% function tgrid = interpolFtoT(fgrid, mask, M, N)'
tgrid = .25 * mask .* ( ...
              fgrid(1:M,1:N) ...
            + fgrid(2:M+1,1:N) ...
            + fgrid(2:M+1,2:N+1) ...
            + fgrid(1:M,2:N+1) ...
            );
'%--------------------------------------------------------------------------
%SHEL SHallow-water numerical modEL
%Copyright (C) 2006,2009,2010. Guillaume Riflet, Instituto Superior Técnico
%da Universidade Técnica de Lisboa.
%--------------------------------------------------------------------------'
```

Right after creating the function, call it from the main function in the computeDiagnostic.m file, [right before the interpolation of U and V to T-cells and the 1D variables computation](http://code.google.com/p/shel/source/browse/trunk/ComputeDiagnostics.m#176):

```
'%%%%Begin new code here%%%%%'
'%dummy1'
computeDummy1_w;
'%%%%End new code%%%%%%%%%%%%'

'%%interpolate from U and V to T cells'
u_t = .5 * ( u(1:M,:) + u(2:M+1,:) );
v_t = .5 * ( v(:,1:N) + v(:,2:N+1) );

Hnoland = H;
Hnoland(find(H < -1)) = 0.;

'%Volume
%volume(l) = dA * sum(sum(Hnoland)); %m3
%Perturbed volume'
volume(l) = dA * sum(sum(eta .* mask)); %m3
```

# Update the property list #
The idea is to add the `dummy1` property in the property list of the visualize\_v4 _four-views panel_
![http://content.screencast.com/users/GRiflet/folders/Jing/media/ef4d9113-81eb-4eb1-b6da-805a3bbe1fbc/propertylist.png](http://content.screencast.com/users/GRiflet/folders/Jing/media/ef4d9113-81eb-4eb1-b6da-805a3bbe1fbc/propertylist.png)

The first thing to do, it to get a sense of where would fit the new variable in the list. The most accessible place is not necessarily the last place. If the property is matricial diagnostic quantity, then the variable should go right before the `global volume` property, at the end of the other diagnostic properties. If the property is a domain integrated, time-evolving property, then it should go at the end of global properties, right after `global Okubo-Weiss` and before `Section X-axis`.

**In this case, the `dummy1` property will be inserted at the end of the diagnostic matricial properties, right before the `global volume` property.**

## visualize\_v4 GUI ##
Call the guide from the matlab prompt, to edit the visualize\_v4 panel:
```
 matlab promtp >> guide visualize_v4
```
Then edit both dropdown menus in the  left and the right-property sub-panels: edit the string and inject the property name in the right place in the list. Then save everything.
![http://content.screencast.com/users/GRiflet/folders/Jing/media/226ec21c-3438-40a8-946b-c17467fecf34/AddPropertyInPropertyList.png](http://content.screencast.com/users/GRiflet/folders/Jing/media/226ec21c-3438-40a8-946b-c17467fecf34/AddPropertyInPropertyList.png)

# Choose the rendering #
For the last part of this tutorial, there will be the need to choose and implement the rendering type of the `dummy1` property.

So far, there are six types of rendering modes:
  1. 3D surface
  1. color map
  1. 2D vector field plot
  1. contour plot
  1. time-series plot
  1. Cross-section plot

**In this case, the `dummy1` property will be rendered using the _3D surface_ mode.**

## plotmodel.m ##

Edit the commented enumerated property list, [here](http://code.google.com/p/shel/source/browse/trunk/plotmodel.m#72):
```
'% 21 - Okubo-Weiss (W)'
'% 22 - Dummy1 (W)'
'% 23 - global volume'
```

Summon the missing global variables, [here](http://code.google.com/p/shel/source/browse/trunk/plotmodel.m#86):
```
    global enstrophy_t;
    global dummy1_t;
    global masknan;
```

Inject the new _case_ in the [switch](http://code.google.com/p/shel/source/browse/trunk/plotmodel.m#145):
```
            case 21 '% Scalar of Okubo-Weiss (1/s2)'
                surf(Hn,x(2:M-1,2:N-1),y(2:M-1,2:N-1), masknan(2:M-1,2:N-1) .* okuboweiss_t(2:M-1,2:N-1));
                '%title(Hn,'Okubo-Weiss');'
                SetXY(Hn);
                SetZ(CLim,'Okubo-Weiss,  \itW\rm (s^{-2})',Hn,az,el);

            case 22 '% Dummy1 (put here the physics dimensions)'
                surf(Hn,x(2:M-1,2:N-1),y(2:M-1,2:N-1), masknan(2:M-1,2:N-1) .* dummy1_t(2:M-1,2:N-1));
                title(Hn,'Dummy1');
                SetXY(Hn);
                SetZ(CLim,'Dummy1,  \itW\rm (physics dimension in latex notation)',Hn,az,el);

            case 23 '%global volume variation(m3)'
                plot(Hn,vtime, volume);
' %               axis(Hn,[min(vtime) max(vtime) min(volume)*.999999 max(volume)*1.000001]);'
                xlabel(Hn,'Time,  \itt\rm (s)', 'FontSize', fontsize);
                ylabel(Hn,'Volume,  \itV\rm (m^3)', 'FontSize', fontsize);
                title(Hn,'Volume deviation', 'FontSize', titfontsize);
                set(Hn, 'FontSize', fontsize);
```

And that's it!

# Now that it's done #

Users can now render the 3D surface of the `dummy1` property. In order to integrate completely the new property in the official SHEL project, please ask the project manager to make you a SHEL developer: `guillaume.riflet at gmail.com`.