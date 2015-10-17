

The templates described here are provided in the google subversion repository and in the software package available to download in the [downloads section](http://code.google.com/p/shel/downloads/list).

---

# Bump description #
<a href='http://www.youtube.com/watch?feature=player_embedded&v=Ibt12ftTAwM' target='_blank'><img src='http://img.youtube.com/vi/Ibt12ftTAwM/0.jpg' width='425' height=344 /></a>

_bump.mat_ file

## Control panel parameters ##
  * uses the _pressure force_ only.
  * zero viscosity.
  * 10 m constant depth.
  * 30 x 30 grid size.
  * 20 m x 20 m spatial increment.
  * 0.1 s time-step, 100 s simulation duration and 3 s output interval.
  * Gaussian bump (refer to the _technical guide_ in the [downloads section](http://code.google.com/p/shel/downloads/list) for more details) with its center at the domain origin (centre). It has a height of 1 cm and a width of 60 m along both x and y-axis.
  * Boundaries: _level_ is set to Sommerfeld, _normal_ is set to _Flather_  and _tangent_ is set to _Dirichelet (IC)_.

## Variants ##
The following variants have the same characteristics than the _bump_ case except when indicated.

### closed bump ###
<a href='http://www.youtube.com/watch?feature=player_embedded&v=pxkpbE-XTaI' target='_blank'><img src='http://img.youtube.com/vi/pxkpbE-XTaI/0.jpg' width='425' height=344 /></a>

  * _closedbump.mat_ file.
  * Boundaries: _u closed_ and _v closed_ are set.

### geostrophic bump ###
<a href='http://www.youtube.com/watch?feature=player_embedded&v=3fK1mvRjR9M' target='_blank'><img src='http://img.youtube.com/vi/3fK1mvRjR9M/0.jpg' width='425' height=344 /></a>

The geostrophic bump uses the Coriolis force also.
  * _geobump.mat_ file.
  * uses the _coriolis force_ also at a _latitude (f)_ of 43 degrees North.
  * spatial-stepping and gaussian width are multiplied by 50 000, making a 100 km x 100 km cell size and a gaussian width of 300 km.
  * time step is of 400 s, total simulation duration is 400000 s and output intervals is 2000.

There is another geostrophic bump configuration for a slightly larger domain:
  * _geobump2.mat_ file.
  * 37 x 37 grid cells.
  * 20 km x 20 km cell size.
  * 60 km x 60 km gaussian width.
  * time-step of 500 s, with a simulation duration of 100000s and an output interval of 1000 s.


---

# Step-bump description #
<a href='http://www.youtube.com/watch?feature=player_embedded&v=lcvYr-M_sk0' target='_blank'><img src='http://img.youtube.com/vi/lcvYr-M_sk0/0.jpg' width='425' height=344 /></a>

_stepbump.mat_ file

## Control panel parameters ##
  * uses the _pressure force_ only.
  * zero viscosity.
  * 10 m depth. Step-bottom with 1 m depth step.
  * Island with 27 km of radius, at a distance from the origin (center of the domain) of 20 km to the South and 80 km to the West.
  * 37 x 37 grid size.
  * 20 km x 20 km spatial increment.
  * 500 s time-step, 100000 s simulation duration and 1000 s output interval.
  * Gaussian bump (refer to the _technical guide_ in the [downloads section](http://code.google.com/p/shel/downloads/list) for more details) with its center at the domain origin (centre). It has a height of 1 cm and a width of 60 km along both x and y-axis.
  * Boundaries: _level_ is set to Sommerfeld, _normal_ is set to _Flather_  and _tangent_ is set to _Dirichelet (IC)_.

## Variants ##
The following variants have the same characteristics than the _step-bump_ case except when indicated.

### Closed step-bump ###
  * _closedstepbump.mat_ file.
  * Boundaries: _u closed_ and _v closed_ are set.

### Geostrophic step-bump ###
<a href='http://www.youtube.com/watch?feature=player_embedded&v=reqIHQgXCq0' target='_blank'><img src='http://img.youtube.com/vi/reqIHQgXCq0/0.jpg' width='425' height=344 /></a>

  * _geostepbump.mat_ file.
  * uses the _coriolis force_ also at a _latitude (f)_ of 43.98 degrees North.


---

# Taylor column description #
_taylor geoIC.mat_ file.

## Control panel parameters ##
  * uses the _pressure force_ and the _Coriolis force_ at a _latitude (f)_ of 43.98 degrees North.
  * 0.5 m/s<sup>2</sup> viscosity.
  * Loaded a bathymetry file from the folder _bathymetries/taylor column flow.mat_, containing the bathymetry variable _`d`_ . The initial velocity field (_`u`_, _`v`_) and water elevation _`eta`_ is over-written by the _Taylor col_ checkbox in the _initial conditions_ sub-panel.
  * 55 x 33 grid size.
  * 10 km x 10 km spatial increment.
  * 300 s time-step, 200000 s simulation duration and 3000 s output interval.
  * _Taylor col_ is checked. This builds the initial geostrophical velocity field and water elevation within the domain.
  * Boundaries: The y-axis is closed (_v-closed_ is checked), making the domain more like a one-way pipeline. The _level_ is set to Sommerfeld, _normal_ is set to _Dirichelet (IC)_  and _tangent_ is set to _Dirichelet (IC)_, meaning that the initial velocity field at the open-boundaries will force continuously the domain throughout the simulation.

## Variants ##
### Null Taylor column ###
<a href='http://www.youtube.com/watch?feature=player_embedded&v=ltTsbICcQEc' target='_blank'><img src='http://img.youtube.com/vi/ltTsbICcQEc/0.jpg' width='425' height=344 /></a>

  * _taylor null IC.mat_ file.
  * _Taylor col_ is left unchecked. This means that the bathymetry file _taylor column flow.mat_ will load, beyond the bathymetry, also the initial velocity field, _`u`_ and _`v`_, and the initial water elevation, _`eta`_. In this configuration, the initial water elevation is null everywhere and the velocity field is null everywhere, except at the open-boundaries, where it has a constant velocity along the x-axis (non-null u-component).


---

# Proposing more templates #
To submit or contribute with new templates, please contact the author: `guillaume.riflet at gmail.com`.