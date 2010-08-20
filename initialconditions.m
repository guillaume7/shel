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

global x;
global y;

global x_u;
global y_u;

global x_v;
global y_v;

global x_w;
global y_w;

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

%T-cells
global Tr;
global eta0;
global eta_old;
global eta;
global eta_new;
global H_old;
global H;
global H_new;
global d;
global mask;
global masknan;
global u_t;
global v_t;
global ke;
global pe;
%Eddy kinetic energy
global eke;
global gradux;
global graduy;
global gradvx;
global gradvy;

global dx;
global dy;

%U-cells
global u_a;
global u_old;
global u;
global u_new;
global mask_u;

%V-cells
global v_a;
global v_old;
global v;
global v_new;
global mask_v;

%Z-cells
global curl_w; %The curl of velocity (v_x - u_y)
global curl_t; % (v_x - u_y) The curl of velocity
global potvorticity_w; %Rossby's potential vorticity
global potvorticity_t;
global strechrate_t; % (u_x - v_y) %Check Arakawa 1966
global shearrate_w;
global shearrate_t; % (v_x + u_y) %Check Arakawa 1966
global sqshearrate_w;
global sqshearrate_t;
global sqstrechrate_t;
global enstrophy_w;
global enstrophy_t; % (0.5 * curl * curl)
global sqstrain_t; % 0.5 * ( shear * shear + strech * strech )
global divergence_t;
global sqdivergence_t;
global okuboweiss_w;
global okuboweiss_t; % (sqstrain - enstrophy) (Check Arakawa1966 and Weiss1981)

%time-dependent global properties
global L;
global volume;
global iKe;
global iPe;
global vtime;
global iVort;
global iEnst;
global iSqStrech;
global iSqShear;
global iSqStrain;
global iOWeiss;
global iMomentumU;
global iMomentumV;
global iEke;

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

    