%--------------------------------------------------------------------------
%     This file is part of SHEL SHallow-water numerical modEL
% 
%     SHEL is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     SHEL is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with SHEL.  If not, see <http://www.gnu.org/licenses/>.
%--------------------------------------------------------------------------
function initialconditions
%function initialconditions

    global time;
    global timetr;
    global dt;
    global frame;

    %GR : clears figure of images
    frame = 1;
    time = 0.;
    timetr = time + 2 * dt;
    updateTime;
    makecoordinates;
    fillfields;

%% makecoordinates
function makecoordinates
%function makecoordinates

global dx;
global dy;
global M;
global N;
global x0;
global y0;

global x; %x-component of the T-grid coordinates. Real matrix.
global y; %y-component of the T-grid coordinates. Real matrix.

global x_u; %x-component of the U-grid coordinates. Real matrix.
global y_u; %y-component of the U-grid coordinates. Real matrix.

global x_v; %x-component of the V-grid coordinates. Real matrix.
global y_v; %y-component of the V-grid coordinates. Real matrix.

global x_w; %x-component of the W-grid coordinates. Real matrix.
global y_w; %y-component of the W-grid coordinates. Real matrix.

global loadbathymetry;
global bathymetryfile;
global varsinfile;

if loadbathymetry
    varsinfile = who('-file',bathymetryfile);
    if ~isempty(varsinfile)
        load(bathymetryfile, 'd');
        siz = size(d);
        M = siz(1);
        N = siz(2);
    end
end

x0 = dx/2;
y0 = dy/2;

%T-cell
x = (1:M)' * (1:N);
y = (1:M)' * (1:N);
for j = 1:N
    x(:,j) = x0 + (x(:,j)/j - 1) * dx;
end
for i = 1:M
    y(i,:) = y0 + (y(i,:)/i - 1) * dy;
end

%U-cell
x_u = (1:M+1)' * (1:N);
y_u = (1:M+1)' * (1:N);
for j = 1:N
    x_u(:,j) = x0 + (x_u(:,j)/j - 1) * dx - dx/2;
end
for i = 1:M
    y_u(i,:) = y0 + (y_u(i,:)/i - 1) * dy;
end

%V-cell
x_v = (1:M)' * (1:N+1);
y_v = (1:M)' * (1:N+1);
for j = 1:N
    x_v(:,j) = x0 + (x_v(:,j)/j - 1) * dx;
end
for i = 1:M
    y_v(i,:) = y0 + (y_v(i,:)/i - 1) * dy - dy/2;
end

%W-cell
x_w = (1:M)' * (1:N);
y_w = (1:M)' * (1:N);
for j = 1:N
    x_w(:,j) = x0 + (x_w(:,j)/j - 1) * dx - dx/2;
end
for i = 1:M
    y_w(i,:) = y0 + (y_w(i,:)/i - 1) * dy - dy/2;
end

%% fillfields
function fillfields
%function fillfields

%         Arakawa C staggered grid (Arakawa 1966)
%
%    --j y
%   |           --- V ---           V ------- V         U -- T -- U
%   i          |         |          |         |         |         |
%   x          U    T    U          T    U    T         |    V    |
%              |         |          |         |         |         |
%               --- V ---           V ------- V         U -- T -- U
%                 T-cell               U-cell              V-cell
%
%       T(eta, H, tr, d, f, mask,   U(u, mask_u         V(v, mask_v
%           x, y, wind and bottom      x_u, y_u)           x_v, y_v)
%            stress)
%    --j y
%   |          T--- U ---T 
%   i          |         | 
%   x          V    W    V 
%              |         | 
%              T--- U ---T 
%                 W-cell   
%
%       W(zeta - curl of v, vorticity)

%Test-cases
global tc_bump;
global tc_geostrophic;
global tc_isla;
global tc_taylor;

%boundary condition
global u_closed;
global v_closed;
%global bc_closed;
global noslip;

%Step bathym
global step;

%bump
global bump_x0;
global bump_y0;
global bump_sx;
global bump_sy;
global bump_d0;

%island
global isla_x0;
global isla_y0;
global isla_R;

%tracer
global tracer;
global TrXo;
global TrYo;
global TrR;

%Grid and bathymetry related
global land;
global M;
global N;
global d0;
global d0_step
global dx;
global dy;

%T-cells
global Tr; % tracer concentration.
global eta0; % initial water elevation.
global eta_old; %water elevation at instant t-dt.
global eta; %water elevation at instant t.
global eta_new; %water elevation at instant t+dt.
global H_old; %water column height at instant t-dt.
global H; %water column height at instant t.
global H_new; %water column height at instant t+dt.
global d; %bathymetry depth. constant in time.
global mask; %land mask.
global masknan; %not-a-number mask.
global u_t; %u component of velocity.
global v_t; %v component of velocity.
 %Diagnostic quantities
global ke; %kinetic energy.
global pe; %potential energy.
global eke; %turbulent kinetic energy.
global gradux; %u gradient along the x-coordinate.
global graduy; %u gradient along the y-coordinate.
global gradvx; %v gradient along the x-coordinate.
global gradvy; %v gradient along the y-coordinate.
global curl_t; % (v_x - u_y) The curl of velocity
global potvorticity_t; % the potential vorticity
global strechrate_t; % strech rate. (u_x - v_y) %Check Arakawa 1966
global shearrate_t; % shear rate. (v_x + u_y) %Check Arakawa 1966
global sqstrechrate_t; %the square of the strech rate.
global sqshearrate_t; %the square of the shear rate.
global enstrophy_t; % enstrophy (0.5 * curl * curl)
global sqstrain_t; % the square of the strain rate. 
                % 0.5 * ( shear * shear + strech * strech )
global divergence_t; %the horizontal divergence.
global sqdivergence_t; %the square of the horizontal divergence.
global okuboweiss_t; % the okubo-weiss scalar (sqstrain - enstrophy) 
                    %(Check Arakawa1966 and Weiss1981)

%U-cells
global u_a; %time averaged u component of velocity.
global u_old; %u component of velocity at instant t-dt.
global u; %u component of velocity at instant t.
global u_new; %u component of velocity at instant t+dt.
global mask_u; %flux mask of the u-component of velocity.

%V-cells
global v_a; %time averaged v component of velocity.
global v_old; %v component of velocity at instant t-dt.
global v; %v component of velocity at instant t.
global v_new; %v component of velocity at instant t+dt.
global mask_v; %flux mask of the v-component of velocity.

%Z-cells
global curl_w; %The curl of velocity (v_x - u_y)
global potvorticity_w; %Rossby's potential vorticity
global shearrate_w; % shear rate. (v_x + u_y) %Check Arakawa 1966
global sqshearrate_w; % the square of the shear rate.
global enstrophy_w; % enstrophy (0.5 * curl * curl)
global okuboweiss_w; % the okubo-weiss scalar (sqstrain - enstrophy) 
                    %(Check Arakawa1966 and Weiss1981)

%time-dependent global properties (domain integrated)
global L;
global volume; %global volume.
global iKe; % kinetic energy.
global iPe; % potential energy.
global vtime; % time vector.
global iVort; % relative vorticity.
global iEnst; % enstrophy.
global iSqStrech; %square of the strech rate.
global iSqShear; %square of the shear rate.
global iSqStrain; %square of the strain rate.
global iOWeiss; %okubo-weiss parameter.
global iMomentumU; %the u-component of momentum. 
global iMomentumV; %the v-component of momentum.
global iEke; % turbulent kinetic energy.

global loadbathymetry;
global bathymetryfile;
global varsinfile;

%Depth field (bathymetry)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
land = -99;

%%Constant depth or read from bathymetry file
if tc_taylor %Submarine mount equals 30% of total mean depth
    d = d0 - cylindermap(zeros(size(d)), 0.3 * d0, bump_sx, bump_sy, ...
                            M * dx * .5 - bump_x0, ...
                            N * dy * .5 - bump_y0);
elseif loadbathymetry && sum(strcmp(varsinfile,'d'))
    load(bathymetryfile, 'd');
else
    d = d0 * ones(M,N);
end

%%changing depth
if step
    d = makestep(d,d0_step);
end

%%Closed boundaries: boundaries are filled with land
if v_closed
    for j = 1:N-1:N;
    for i = 1:M
        d(i,j) = land;
    end
    end
end

if u_closed    
    for i = 1:M-1:M;
    for j = 1:N
        d(i,j) = land;
    end
    end
end

%%Inner islands (random)
%%%arguments: depth, io, jo, r
if tc_isla
    d = makeisland(d, M * dx * .5 - isla_x0, N * dy * .5 - isla_y0, isla_R);
end

%WATERMASK matrices%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%T-cell
%%all-water
mask = ones(M,N);

%%eventual islands and closed domain
for j=1:N
for i=1:M
    if d(i,j) < 0.
        mask(i,j) = 0;
    end
end
end

%U-Cell and V-cell
mask_u = ones(M+1,N);
mask_v = ones(M,N+1);
for j = 1:N
for i = 1:M
    if mask( i, j) == 0
        mask_u(  i,  j) = mask(i,j);
        mask_u(i+1,  j) = mask(i,j);
        mask_v(  i,  j) = mask(i,j);
        mask_v(  i,j+1) = mask(i,j);
    end
end
end

if noslip
    %U-cells
    for j = 2:N-1
    for i = 2:M
        mask_u(i,j) = mask(i  ,j+1) * mask(i  ,j-1) * mask(i-1,j+1) * mask(i-1,j-1);
    end
    end
    %V-Cells
    for j = 2:N
    for i = 2:M-1
        mask_v(i,j) = mask(i+1,  j) * mask(i+1,j-1) * mask(i-1,  j) * mask(i-1,j-1);
    end
    end
end

%Water level fields%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (tc_bump || tc_geostrophic)
    eta_old = bumpmap(bump_d0, bump_sx, bump_sy, M * dx * .5 - bump_x0, ...
                                                N * dy * .5 - bump_y0);
else
    if loadbathymetry && sum(strcmp(varsinfile,'eta'))
        load(bathymetryfile, 'eta');
        eta_old = eta;
    else
        eta_old = eta0 * ones(M,N);
    end
end
eta = eta_old;
eta_new = eta;

%Bottom-to-surface depth (always positive below the water surface)
%H = eta + d.
H_old = eta_old + d;
H =H_old;
H_new = H_old;

if tracer
    %Unitary concentration...
    Tr = zeros(M,N);
    Tr = cylindermap(Tr, 1., TrR, TrR, M * dx * .5 - TrXo, ...
                                    N * dy * .5 - TrYo);
end

%U velocity fields%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if tc_geostrophic
    u_old = ugeo(eta);    
elseif tc_taylor % 3 cm/s u-flow
    u_old = zeros(size(mask_u));
    %%uncomment one of the following lines
    u_old = 0.03 * ones(size(mask_u));
    %u_old(1:M:M+1,:) = 0.03 * ones(size(mask_u(1:M:M+1,:)));
else
    if loadbathymetry && sum(strcmp(varsinfile,'u'))
        load(bathymetryfile, 'u');
        u_old = u;
    else
        u_old = zeros(M+1,N);
    end
end
u = u_old;
u_new = u;
u_a = u;

%V velocity fields%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if tc_geostrophic
    v_old = vgeo(eta);    
else
    if loadbathymetry && sum(strcmp(varsinfile,'v'))
        load(bathymetryfile, 'v');
        v_old = v;
    else
        v_old = zeros(M,N+1);
    end
end
v = v_old;
v_new = v;
v_a = v;

%call this function when you need 
%to manually initialize the level
%in the Taylor column flow.
if tc_taylor
    makeGeostLevelFromU
end

%Z field%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
curl_w = zeros(M+1,N+1);
potvorticity_w = zeros(M+1,N+1);
enstrophy_w = zeros(M+1,N+1);
shearrate_w = zeros(M+1,N+1);
sqshearrate_w = zeros(M+1,N+1);
okuboweiss_w = zeros(M+1,N+1);

%Visualization matrices%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
masknan = ones(M,N);
for j=1:N
for i=1:M
        if mask(i,j) == 0
            masknan(i,j) = nan;
        else
            masknan(i,j) = 1;
        end
end
end

u_t = zeros(M,N);
v_t = zeros(M,N);

ke = zeros(M,N);
pe = zeros(M,N);

curl_t = zeros(M,N);
potvorticity_t = zeros(M,N);
strechrate_t = zeros(M,N);
shearrate_t = zeros(M,N);
enstrophy_t = zeros(M,N);
sqstrain_t = zeros(M,N);
okuboweiss_t = zeros(M,N);
sqshearrate_t = zeros(M,N);
sqstrechrate_t = zeros(M,N);
divergence_t = zeros(M,N);
sqdivergence_t = zeros(M,N);

%Eddy kinetic energy
gradux=zeros(M,N);
graduy=zeros(M,N);
gradvx=zeros(M,N);
gradvy=zeros(M,N);
eke=zeros(M,N);

%One-dimensional diagnostic variables
%Energy initialization
vtime = nan(1,L);
iKe = nan(1,L);
iPe = nan(1,L);
iEke = nan(1,L);
volume = nan(1,L);
iVort = nan(1,L);
iEnst = nan(1,L);
iSqShear = nan(1,L);
iSqStrech = nan(1,L);
iSqStrain = nan(1,L);
iOWeiss = nan(1,L);
iMomentumU = nan(1,L);
iMomentumV = nan(1,L);

ComputeDiagnostics(1);

%cylinder field
function map_L = cylindermap(map_L, l, sx, sy, x_L, y_L)

    global x;
    global y;
       
    ind = find ( (x - x_L).^2 / sx^2 + (y - y_L).^2 / sy^2 < 1. );
    map_L(ind) = l;
    
%bump test-case field
function map_L = bumpmap(l, sx, sy, x_L, y_L)

    global M;
    global N;
    global x;
    global y;
       
    map_L = zeros(M,N);   
    
    for i = 1:M
    for j = 1:N
        map_L(i,j) = l * exp( - ( ( x(i,j) - x_L )^2 / sx^2 + ( y(i,j) - y_L )^2 / sy^2 ) );
    end
    end
        
function depth = makeisland(depth, x_L, y_L, r)
%function depth = makeisland(depth, x_L, y_L, r)

    global M;
    global N;
    global x;
    global y;
    global land;
    
    for i = 1:M
    for j = 1:N
        if r^2 > ( x(i,j) - x_L )^2 + ( y(i,j) - y_L )^2
            depth(i,j) = land;
        end
    end
    end

function ugeo_L = ugeo(eta_L)
%function ugeo_L = ugeo(eta_L)

    global M;
    global N;
    global g;
    global f;
    
    j_L = floor(N/2)+1;
    sig_L = 1.;
    ug_L = zeros(M,N);
    
    for i=1:M
    for j=1:N
        ug_L(i,j) = g / f * ( j - j_L ) / sig_L^2 * eta_L(i,j);
    end
    end
    ugeo_L = InterpolTU(ug_L);

function vgeo_L = vgeo(eta_L)
%function vgeo_L = vgeo(eta_L)
    global M;
    global N;
    global g;
    global f;
    i_L = floor(M/2)+1;
    sig_L = 1.;
    vg_L = zeros(M,N);
    
    for i=1:M
    for j=1:N
        vg_L(i,j) = - g / f * ( i - i_L ) / sig_L^2 * eta_L(i,j);
    end
    end
    vgeo_L = InterpolTV(vg_L);

function bath = makestep(bath,depth)
%function bath = makestep(bath,depth)
%Creates a stepped bathymetry
    sz = size(bath);
    M = sz(1);
    N = sz(2);
    j0 = N/5;
    
    for i = 1:M
    for j = 1:N
        if j < 0.9 * i + j0
            bath(i,j) = depth;
        end
    end
    end

function makeGeostLevelFromU

    global eta_old
    global eta
    global eta_new
    global d
    global H_old
    global H
    global H_new
    global f
    global g
    global mask
    global mask_u
    global dy
    global u
    global M
    global N

    for j=2:N
    for i=1:M
        eta(i,j) = eta(i,j-1) * mask(i,j-1) - dy * f/g * (...
            mask_u(i,j) * u(i,j) + mask_u(i+1,j) * u(i+1,j) ...
            + mask_u(i,j-1) * u(i,j-1) + mask_u(i+1,j-1) * u(i+1,j-1) ...
            ) / ( ...
            mask_u(i,j) + mask_u(i+1,j) + mask_u(i,j-1) + mask_u(i+1,j-1) ...
            );
        eta(i,j) = mask(i,j) * eta(i,j);
    end
    end

    eta_old = eta;
    eta_new = eta;

    H_old = eta_old + d;
    H =H_old;
    H_new = H_old;


%--------------------------------------------------------------------------
%SHEL SHallow-water numerical modEL
%Copyright (C) 2006,2009,2010. Guillaume Riflet, Instituto Superior Técnico
%da Universidade Técnica de Lisboa.
%--------------------------------------------------------------------------
