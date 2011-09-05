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

function s = model_handles2
%function s = model_handles2
%
%Call it this way:
% s = model_handles;
% y = s.model_compute_input(); % --> the parenthesis ARE REQUIRED!
% y = s.model_compute_output( handles );

s.model_compute_input = @ComputeModelInput;
s.model_compute_output = @ComputeModelOutput;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Here starts the model Inputs part
%%Initial conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ComputeModelInput
%function ComputeModelInputs

    global time;
    global timetr;
    global dt;
    global frame;

    %GR : clears figure of images
    frame = 1;
    time = 0.;
    timetr = time + 2 * dt;
    ComputeTime;
    makecoordinates;
    fillfields;

function ComputeTime
%function updateTime

    global dt;
    global duration;
    global outputdt;
    global L;
    global outputL;

    L = getL(duration,dt);
    outputL = getL(outputdt, dt);

function L_L = getL (duration_L, dt_L)
%function L_L = getL (duration_L, dt_L)

    L_L = floor(duration_L / dt_L);
    
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Here starts the model Outputs part
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ComputeModelOutput( handles )
%function ComputeModelOutputs( handles ) [old ComputeModel]

%Time parameters and output parameters
global L;
global outputL;
%global dt;
global time;
%global timetr;
global frame;
global printG;
global film;
global mov;
global myfullfile;
global mydir;
global optprint;

%tracer parameters
global tracer;
%global Tr;
%global K;

s = GUI_common_handles;

%Movie initialization
clear mex;
tic;
ol = 0;
movie = false;

%Define the type of output
switch optprint
    case {11, 12, 13, 14, 15}
        movie = true;
end

% Write output of initial condition to disk.
if printG    
    if movie
        %Open a new avi file.
        filename = sprintf('%s.avi', [mydir,'/',myfullfile]);
        mov = avifile(filename,'fps', 10,'quality',100,'compression','None');
    end
    s.GUI_printit(frame, handles, optprint);
    frame = frame + 1;
end 

% Main loop in time
ltr = 0;
for l = 1 : L           
    
    %Takes care of the numerical scheme
    %for momentum and for the tracer
    %This line can easily be replaced
    %with a switch, choosing one kernel
    %among many (mex file kernel, other
    %kernels ...).
    kernel_SHEL(l);
    
% Take care of plotting the outputs and write them to disk.
    if ol == outputL
        
        ol=0;
               
        if film
            s.GUI_plotmodel(handles);
            pause(0.1); %s
        end
        if printG
            s.GUI_printit(frame, handles, optprint);
            frame = frame + 1;
        end

        %Display time
        strin = sprintf('time = %0.3g s', time);
        set(handles.timetext,'String',strin);
        pause(0.01); %s pause (required if I want to update the clock).
        
    end
    ol=ol+1;
    
%End of main loop in time    
end

% Take care of writing the final output on disk
if printG
if movie
    mov = close(mov);
end
end

if ~film
    s.GUI_plotmodel(handles);
end

% Display simulation statistics and performance
cput = toc;
disp(sprintf('Elapsed time: %0.5f s\nTime per iteration: %0.5f s\n', cput, cput/L));
if tracer
    disp(sprintf('Made %0.3d tracer iterations on a total of %0.3d iterations', ltr, l));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Here starts the Kernel part
%%the numerical scheme
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function kernel_SHEL(l)
%function kernel_SHEL (old ComputeModel)

%Explicit leapfrog scheme
%combined with Asselin-Roberts filter
%(the filter supresses the computational mode
%and prevents decoupling between odd/even timesteps)
%
%         Arakawa C staggered grid (Arakawa 1966)
%
%    --j y
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
%        (deformation rate = curl(-u,v,0)
%%%

%Time parameters and output parameters
global dt;
global time;
global timetr;

%tracer parameters
global tracer;
global Tr;
global K;
   
    %Update the time
    time = time + dt;
    
    %Compute the time scheme of the momentum+continuity model
    ComputeLeapfrog;

    %Compute all the tracers
    if tracer
        if timetr <= time
            [Tr, dttr] = ComputeTracer_FT(Tr, K);
            timetr = time + dttr;
            ltr = ltr + 1;
        end
    end
    
    %Numerics to compute the diagnostic variables
    ComputeDiagnostics(l);

function ComputeLeapfrog
%Shallow-water equations solver
%
global gama;
global dx;
global dy;
global dt;
 
%T-cells
global eta_old;
global eta;
global eta_new;
global H_old;
global H;
global H_new;
global d;
global mask;
global dtmask;
 
%U-cells
global u_old;
global u;
global u_new;
global mask_u;
global dtmask_u;
 
%V-cells
global v_old; 
global v;
global v_new;
global mask_v;
global dtmask_v;

%% Solve the spatial schemes in the interior of the domain
    % For the water elevation
    RHSeta = ComputeContinuity(H, eta, u, v, mask);
    eta_new = eta_old + 2 * dt * RHSeta .* mask; 
    H_new = eta_new + d;
    ind = find(H_new < 0);
    eta_new(ind) = -d;
    dtmask = (eta_new - eta_old) * .5 ./ RHSeta .* mask;
    %Must build the dtmask_u and dtmask_v in order to preserve momentum
    %too.
    dtmask_u = buildmask_U(dtmask,mask_u);
    dtmask_v = buildmask_V(dtmask',mask_v')';
    
    RHSeta = ComputeContinuity(H, eta, u, v, dtmask);
    
    % For the u and v components of velocity
    % % The old 10%-faster way (but a lot more error-prone)
    %     [RHSu, RHSv] = ComputeUV2D_old;
    % Switch v <-> u; v_old <-> u_old; dy <-> dx and -1. <-> 1. (coriolis)
    % Transpose the matrices.
    RHSu = ComputeSpaceU_CS( H, H_old, eta, ...
                        u, u_old, v, v_old, ...
                        dtmask_u, dx, dy, 1.);
    RHSv = ComputeSpaceU_CS( H', H_old', eta', ...
                        v', v_old', u', u_old', ...
                        dtmask_v', dy, dx, -1.)';
 
%% LEAPFROG TIME-SCHEME
    % Water level time scheme
    eta_new = eta_old + 2 * RHSeta .* mask; 
    H_new = eta_new + d;
    %u and v time schemes inside the domain. 
    %Leapfrog means that 2xdt is passed as teh time step.
    v_new = ComputeTimeU_FT(v_old', v_new', RHSv', H_old', H_new', mask_v', 2.)';
    u_new = ComputeTimeU_FT(u_old, u_new, RHSu, H_old, H_new, mask_u, 2.);

%% time scheme at the northern/southern and eastern/southern OB (2 x d
%% <- LeapFrog)

    %Southern and northern boundaries (along V -> transpose inputs)
    %(transposed inputs -> transposed outputs).
    [eta_new, H_new, v_new, u_new] = ComputeOB_U( ...
            eta_new', eta_old', H_new', H_old', d',...
            v_new', mask_v', ...
            u_new', u_old', mask_u', ...
            2*dt, dy ...
     );
    %Eastern and western boundaries (along U)
    %The previous outputs were transposed, so we transpose
    %them back as inputs.
    [eta_new, H_new, u_new, v_new] = ComputeOB_U( ...
            eta_new', eta_old, H_new', H_old, d,...
            u_new', mask_u, ...
            v_new', v_old, mask_v, ...
            2*dt, dx ...
     );
     
%% Asselin-Robert filter 
    eta = eta + gama * ( eta_old - 2 * eta + eta_new ); 
    u = u + gama * ( u_old - 2 * u + u_new ); 
    v = v + gama * ( v_old - 2 * v + v_new ); 

%% Update matrices
    eta_old = eta;
    eta = eta_new;
     
    H_old = H;
    H = H_new;
 
    u_old = u;
    u = u_new;
     
    v_old = v;
    v = v_new;

function RHSeta = ComputeContinuity(H, eta, u, v, mask)
%function ComputeContinuity
%Shallow-water equations solver

% global H; 
% global eta; 
% global u; 
% global v; 
% global mask; 
% global RHSeta;
global dx; 
global dy; 

[M,N] = size(eta);
RHSeta = zeros(M,N); 

%Compute along U then along V
RHSeta(2:M-1,2:N-1) = - ( ...
    ComputeContinuityU(H,u,mask,dx) + ComputeContinuityU(H',v',mask',dy)' ...
);

function RHS = ComputeContinuityU(H,u,mask,dx)
[M,N] = size(H);
RHS = ( ... 
            mask(3:M,2:N-1) .* .5 ...
            .* ( H(3:M,2:N-1) + H(2:M-1,2:N-1) ) .* u(3:M,2:N-1) ... %i+1, i, i+1 
            - mask(1:M-2,2:N-1) .* .5 ...
            .* ( H(2:M-1,2:N-1) + H(1:M-2,2:N-1) ) .* u(2:M-1,2:N-1) ... %i, i-1, i 
        )/dx;

function RHSu = ComputeSpaceU_CS( ...
                    H, H_old, eta, ...
                    u, u_old, v, v_old, ...
                    mask_u, dx, dy, ...
                    signcoriolis ...
                )
%% function newu = ComputeSpaceU_CS
%Shallow-water equations solver for the U component
%of the flow Centred in Space (CS).
%
% Centred differences (CD) or Centred in Space (CS).
%
%NOTE: Any viscous and frictional term
%must be evaluated in tindex(i-1,M)e l-1.
%
% switch v <-> u; j <-> i and dy <-> dx relative to ComputeU2D
% switch v_old <-> u_old
% take special care for the coriolis force term

%% global variables
global g; 
global f; 
global Nu; 
global wind; 
global bottom; 
global pressure; 
global coriolis; 
global wall;

global uwind; 
global vwind; 
global rho0;
global rho_air;
global lb; 
global karman; 

wall = false; 

%% U spatial scheme
[M, N] = size(eta);
RHSu = zeros(size(u));
H_u = zeros(size(u));
H_old_u = zeros(size(u));
v_old_u = zeros(size(u));
v_u = zeros(size(u));
H_u(2:M,2:N-1) = twoaverage_u(H,M,N);
H_old_u(2:M,2:N-1) = twoaverage_u(H_old,M,N);
v_old_u(2:M,2:N-1) = fouraverage_u(v_old(:,2:N),M, N-2); 
v_u(2:M,2:N-1) = fouraverage_u(v(:,2:N),M, N-2); 

RHSu(2:M,2:N-1) = mask_u(2:M,2:N-1) .* ( ... 
    ...%advection along U ... 
        - 0.2500 / dx * ( ... %i-1, i 
        mask_u(3:M+1,2:N-1) .* ( u(3:M+1,2:N-1) + u(2:M,2:N-1) ).^2 .* H(2:M,2:N-1) ... %i+1, i, i+1, i, i 
      - mask_u(1:M-1,2:N-1) .* ( u(2:M,2:N-1) + u(1:M-1,2:N-1) ).^2 .* H(1:M-1,2:N-1) ... %i, i-1, i, i-1, i-1 
     )  ... 
    ...%advection along V 
     - 0.0625 / dy * ( ... %i-1, i 
        ( ... 
            mask_u(2:M,3:N) .* ( H(2:M,3:N) + H(2:M,2:N-1) + H(1:M-1,3:N) + H(1:M-1,2:N-1) ) ... % j+1 ; j ; i-1,j+1 ; i-1 
          .* ( v(2:M,3:N) + v(1:M-1,3:N) ) .* ( u(2:M,3:N) + u(2:M,2:N-1) ) ...% j+1 ; i-1,j+1 ; j+1 ; i 
          - mask_u(2:M,1:N-2) .* ( H(2:M,2:N-1) + H(2:M,1:N-2) + H(1:M-1,2:N-1) + H(1:M-1,1:N-2) ) ... % i ; j-1 ; i-1 ; i-1,j-1 
          .* ( v(2:M,2:N-1) + v(1:M-1,2:N-1) ) .* ( u(2:M,2:N-1) + u(2:M,1:N-2) ) ...% i ; i-1,j-1 ; i ; j-1 
        ) ... 
     ) ... 
    ...%Coriolis Force (U -> +; V -> -)
    + coriolis * ( ... 
        signcoriolis * f * v_u(2:M,2:N-1) .* H_u(2:M,2:N-1) ... % i ; i-1 . 
    )... 
    ...%Pressure gradient Force 
    + pressure * ( ... 
       - g * H_u(2:M,2:N-1) .* ( eta(2:M,2:N-1) - eta(1:M-1,2:N-1) ) / dx ...  
                               ... % i ; i-1 ;                 i ; i-1 
    )... 
    ...%Viscosity along U 
      + ( ... % i ; i-1 
          + Nu * mask_u(3:M+1,2:N-1) .* H_old(2:M,2:N-1) .* ( u_old(3:M+1,2:N-1) - u_old(2:M,2:N-1) ) / dx ... % i ; i+1 ; i ; i 
          - Nu * mask_u(1:M-1,2:N-1) .* H_old(1:M-1,2:N-1) .* ( u_old(2:M,2:N-1) - u_old(1:M-1,2:N-1) ) / dx ... % i-1 ; i ; i-1; i-1 
      ) / dx ... % i ; i-1 
    ...%Viscosity along V 
      + ( ... % i ; i-1 ; i ; i-1 . 
          + Nu / dy * mask_u(2:M,3:N) .* H_old(2:M,2:N-1) .* ( u_old(2:M,3:N) - u_old(2:M,2:N-1) ) ... % (j;j+1;i-1,j+1;i-1);j+1 ; j 
          - Nu / dy * mask_u(2:M,1:N-2) .* H_old(2:M,1:N-2) .* ( u_old(2:M,2:N-1) - u_old(2:M,1:N-2) ) ... % (j-1;j;i-1;i-1,j-1);j ; j-1 
      ) / dy ...
);
 
%Add bottom stress. Look out, so that the bottom stress doesn't
%reverts the direction of momentum ... This wouldn't be very realistic.
if bottom || wall 
    stress = zeros(size(RHSu(2:M,2:N-1))); 
    if bottom 
        stress = bottomstress( u_old(2:M,2:N-1), v_old_u(2:M,2:N-1), ...
                                H_old_u(2:M,2:N-1), ...
                                lb, karman ); 
    end 
    if wall 
        stress = stress + bottomstress(u_old(2:M,2:N-1), 0, dy, lb, karman) ...
            .* ( ... 
                    (1 - mask(1:M-1,3:N)) .* H_old(1:M-1,3:N) + ... 
                    (1 - mask(2:M,3:N)) .* H_old(2:M,3:N) + ... 
                    (1 - mask(1:M-1,1:N-2)) .* H_old(1:M-1,1:N-2) + ... 
                    (1 - mask(2:M,1:N-2)) .* H_old(2:M,1:N-2) ... 
            ) * 0.5 / dy; 
    end 
    ind = find( (RHSu(2:M,2:N-1) - mask_u(2:M,2:N-1) .* stress) ./ RHSu(2:M,2:N-1) > 0 );
    RHSu(ind) = RHSu(ind) - mask_u(ind) .* stress(ind); 
    ind = find( (RHSu(2:M,2:N-1) - mask_u(2:M,2:N-1) .* stress) ./ RHSu(2:M,2:N-1) <= 0 );
    RHSu(ind) = 0.; 
end
    
%Add wind stress.
if wind 
    RHSu(2:M,2:N-1) = RHSu(2:M,2:N-1) ...
            + mask_u(2:M,2:N-1) .* windstress(...
                uwind * ones(size(u(2:M,2:N-1))), ...
                vwind * ones(size(u(2:M,2:N-1))), ...
                rho0, rho_air ...
            );
end
    
function u_new = ComputeTimeU_FT(u_old, u_new, RHSu, H_old, H_new, mask_u, dt)
%function u_new = ComputeTimeU_FT(u_old, u_new, RHSu, H_old, H_new, mask_u, dt)        
%
%Numerical solver with the Forward in Time (FT) Euler time stepping method.
%
    [M,N] = size(H_old);
    
    u_new(2:M,2:N-1) = mask_u(2:M,2:N-1) .* ( ... 
                u_old(2:M,2:N-1) .* ( H_old(2:M,2:N-1) + H_old(1:M-1,2:N-1) ) ... 
                + 2 * dt * RHSu(2:M,2:N-1) ... 
              ) ... 
              ./ ( H_new(2:M,2:N-1) + H_new(1:M-1,2:N-1) );
          
function [eta_new, H_new, u_new, v_new] = ComputeOB_U( ...
    eta_new, eta_old, H_new, H_old, d,...
    u_new, mask_u, ...
    v_new, v_old, mask_v, ...
    dt, dx ...
)
%function computeOB
%Solve the north, south, east and west boundary condition
%for the eta, u and v variables.

global g; 
  
%BC conditions 
global flather;
global neumann;
global radiatetan; 
global radiatelevel;

flather = 0;
neumann = 0;
radiatetan = 0; 
radiatelevel = 1;

[M,N] = size(d);

%% Eastern (1,:) and Western (M,:) boundaries    
if radiatelevel 
    %Waterlevel: GW radiation 
    eta_new(1:M-1:M,:) = mask_u(1:M:M+1,:) .* ( ...
        eta_old(1:M-1:M,:)- dt / dx .* sqrt( g * H_old(1:M-1:M,:) ) ... 
        .* ( eta_old(1:M-1:M,:) - eta_old(2:M-3:M-1,:) ) ...
    ); 
    %no mask below, cuz H must yield -99 where land is.W
    H_new(1:M-1:M,:) = eta_new(1:M-1:M,:) + d(1:M-1:M,:);
end 
if flather 
    %Normal velocity: Flather 
    u_new(1:M:M+1,:) = [-1 * ones(1,N); ones(1,N)] .*  mask_u(1:M:M+1,:) ...
        .* sqrt( g ./ H_new(1:M-1:M,:) ) .* eta_new(1:M-1:M,:); 
elseif neumann
    %Normal velocity: Null-gradient 
    u_new(1:M:M+1,1:N) = mask_u(1:M:M+1,1:N) .* u_new(2:M-2:M,1:N); 
end 
if radiatetan 
    %Tangential velocity: GW radiation 
    v_new(1:M-1:M,2:N) = mask_v(1:M-1:M,2:N) .* ... 
        ( ... 
            v_old(1:M-1:M,2:N) .* ( H_old(1:M-1:M,2:N) + H_old(1:M-1:M,1:N-1) ) ... 
            - dt / dx .* sqrt( g * .5 * ( H_old(1:M-1:M,2:N) + H_old(1:M-1:M,1:N-1) ) ) ... 
            .* ( v_old(1:M-1:M,2:N) - v_old( 2:M-3:M-1,2:N) ) ... 
        ) ... 
        ./ ( H_new(1:M-1:M,2:N) + H_new(1:M-1:M,1:N-1) ); 
end


function [Tr, dttr] = ComputeTracer_FT(Tr, K_L)
%Advection-diffusion tracer equation solver
%Explicit forward in time and upwind
%
global dx;
global dy;
global dt;
 
%T-cells
global H;
global H_new;
global mask;
 
%U-cells
global u_a;
global mask_u;

%V-cells
global v_a;
global mask_v;

%% Solve the tracer spatial scheme in the whole of the domain
    RHSt = mask .* ( ...
        ComputeSpaceT_UP_U(H,u_a,mask_u, Tr, K_L, dx) ... %Along U
       + ComputeSpaceT_UP_U(H',v_a',mask_v',Tr',K_L,dy)' ... %then along V
    );

%% FIND THE OPTIMAL TRACER DT AS A MULTIPLE OF MOMENTUM DT
% within the reasonable limit of maximum 100 * dt.
    minRHSt  = min(min(RHSt));
    ind = find(RHSt == minRHSt);
    dttr = min( floor( abs(- H(ind) .* Tr(ind) ./ minRHSt) / dt ), 1000) ...
        * dt;

%% TRACER FORWARD-IN-TIME SCHEME
    Tr_new = mask .* (H .* Tr + dttr * RHSt) ./ H_new;

%% Update matrices
    Tr = Tr_new;

function RHSt = ComputeSpaceT_UP_U(H, u, mask_u, Tr, K_L, dx_L)
%function RHSt = ComputeSpaceT_UP_U(H, u, mask_u, Tr, K, dx)
%
%UPWIND
%
    global nullgrad;
    nullgrad = false;
    
    [M,N]=size(H);
    Hm = zeros(M+1,N);
    TrUm = zeros(M+1,N);
    TrUp = zeros(M+1,N);
    TrDif = zeros(M+1,N);

    %Computing quantities for the interior of the domain
    TrUm(1:M,:) = .5 * (abs(u(1:M,:)) - u(1:M,:)) .* Tr;
    TrUp(2:M+1,:) = .5 * (abs(u(2:M+1,:)) + u(2:M+1,:)) .* Tr;
    TrDif(2:M,:) = Tr(1:M-1,:) - Tr(2:M,:);
    Hm(2:M,:) = .5 * ( H(1:M-1,:) + H(2:M,:) );
    
    %Computing the quantities for the boundaries
    if nullgrad %Considering null-gradient:
        TrUm(M+1,:) = .5 * (abs(u(M+1,:)) - u(M+1,:)) .* Tr(M,:);
        TrUp(1,:) = .5 * (abs(u(1,:)) + u(1,:)) .* Tr(1,:);
        TrDif(1:M:M+1,:) = 0.;
    else %else null value is considered
        TrUm(M+1,:) = .5 * (abs(u(M+1,:)) - u(M+1,:)) * 0.;
        TrUp(1,:) = .5 * (abs(u(1,:)) + u(1,:)) * 0.;
        TrDif(1,:) = - Tr(1,:);
        TrDif(M+1,:) = Tr(M,:);
    end
    
    %This algo considers that the bathymetry has always a null gradient
    %at the boundary:
    Hm(1:M:M+1,:) = H(1:M-1:M,:);

    %Computing the fluxes at each face of the M+1 faces
    FFluxU = ComputeFaceFluxT_UP(Hm, TrUp, TrUm, TrDif, dx_L, K_L);
    
    %The T-cell increment is the flux balance.
    RHSt = mask_u(1:M,:) .* FFluxU(1:M,:) ...
         - mask_u(2:M+1,:) .* FFluxU(2:M+1,:);

function RHSt = ComputeFaceFluxT_UP(Hm, TrUp, TrUm, TrDif, dx_L, K_L)
%function RHSt = ComputeFaceFluxT_UP(Hm, TrUp, TrUm, TrDif, dx_L, K_L)
%Pure upwind scheme.
    RHSt = Hm .* ( ...
            TrUp - TrUm ...
            + K_L / dx_L * TrDif ...
    ) / dx_L;
        
function av = fouraverage_u(v, M, N) 
%% function av = fouraverage_u(v, mask, M, N)

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

function tb_L = bottomstress( u_L, v_L, H_L, lb, karman) 
% Bottom stress        
%function tb_L = bottomstress( u_L, v_L, H_L ) 
% 
%H_L is the total depth height 
     
    cd_L = karman / log( ( 0.5 .* H_L + lb ) / lb ) ; 
    % Checked-out Blumberg and Mellor 1987 
    % Double-checked with Pietrzak 2002 
    cd_L = max(cd_L .* cd_L, 0.0025); 
    tb_L = cd_L .* sqrt ( u_L .* u_L + v_L .* v_L ) .* u_L; 

function u_L = windstress(u_L, v_L, rho0, rho_air) 
%% Wind stress            
%function tv_L = windstress(v_L) 

    vmod_L = sqrt( u_L .* u_L + v_L .* v_L );

    if vmod_L < 6.
        Cd_L = -.26097 .* vmod_L +  2.4316;
    else
        Cd_L = .064986 .* vmod_L + .44053;
    end

    Cd_L = Cd_L * .001;

    u_L = Cd_L .* rho_air / rho0 .* vmod_L .* u_L;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Here starts the diagnostics part
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function diag = ComputeDiagnostics(l)
%function ComputeDiagnostics
%
%Computes global:
% volume (m3)
% kinetic energy (J)
% potential energy (J)
% total energy (J)
% vorticity (m2/s)
% enstrophy (J)
% squared normal strain (J)
% squared shear strain (J)
% squared strain (J)
% Okubo-Weiss
% momentum (m^3 * m / s)a
%
%Computes local:
% kinetic energy (J)
% potential energy (J)
% total energy (J)
% curl (1/s)
% divergence (1/s)
% stretch (normal strain) (1/s)
% shear (shear strain) (1/s)
% squared strain (1/s)
% enstrophy (1/s2)
% OkuboWeiss (1/s2)
%

global u_a;
global v_a;
global u;
global v;
global u_t;
global v_t;
global H;
global eta;
global M;
global N;
global dx;
global dy;
global dt;
global rho0;
global g;
global K;
global mask;
global ke;
global pe;
global strechrate_t;
global shearrate_w;
global divergence_t;
global sqstrain_t;
global sqstrechrate_t;
global sqshearrate_w;
global curl_w;
global enstrophy_t;
global enstrophy_w;
global okuboweiss_t;

%Eddy kinetic energy
global gradux;
global graduy;
global gradvx;
global gradvy;
global eke;

%time-dependent global properties
global time;
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

dA = dx * dy; %m2

%Computes the average velocity field
u_a = ((l - 1) * u_a + u ) / l;
v_a = ((l - 1) * v_a + v ) / l;

%Eddy kinetic energy
gradux = (u(2:M+1,:) - u(1:M,:)) / dx;
gradvy = (v(:,2:N+1) - v(:,1:N)) / dy;
graduy(2:M,2:N-2) = (u(2:M,3:N-1) - u(2:M,2:N-2)) / dy;
gradvx(2:M-2,2:N) = (v(3:M-1,2:N) - v(2:M-2,2:N)) / dx;

%Compute diagnostic fields
%Curl
computeCurl_w; % 1 / s

%Strain rate
strechrate_t = mask .* ((u(2:M+1,:) - u(1:M,:)) * dy ...
            - (v(:,2:N+1) - v(:,1:N)) * dx) / dA; % 1 / s
computeShearRate_w % 1 / s

%Growth rate
divergence_t = mask .* ((u(2:M+1,1:N) - u(1:M,1:N)) * dy ...
                + (v(1:M,2:N+1) - v(1:M,1:N)) * dx) / dA; % 1/ s

%Potential vorticity
computePotentialVorticity_w

%This is the so-called scalar of Arakawa-Okubo-Weiss.
%Computes the local enstrophy, squared strain rate and Okubo-Weiss
%parameter %Look for Weiss
%La Jolla Institute Report from 1981.
%Physica D: Nonlinear 1991.
%Look for Okubo-Weiss parameter/function
%Look for Arakawa 1966 paper.
enstrophy_w = .5 * ( curl_w.^2); %1 / s2
enstrophy_t = interpolFtoT(enstrophy_w, mask, M, N);

sqshearrate_w = .5 * shearrate_w.^2;
sqshearrate_t = interpolFtoT(sqshearrate_w, mask, M, N);

sqstrechrate_t = .5 * strechrate_t.^2 .* mask;

sqstrain_t = sqstrechrate_t + sqshearrate_t;

sqdivergence_t = .5 * divergence_t.^2 .* mask;

okuboweiss_w = (sqshearrate_w - enstrophy_w);
okuboweiss_t = interpolFtoT(okuboweiss_w, mask, M, N);
okuboweiss_t = okuboweiss_t + sqstrechrate_t .* mask;
%versão do Aires do escalar de Okubo-Weiss
okuboweiss_t = okuboweiss_t  - sqdivergence_t;

%%interpolate from U and V to T cells
u_t = .5 * ( u(1:M,:) + u(2:M+1,:) );
v_t = .5 * ( v(:,1:N) + v(:,2:N+1) );

Hnoland = H;
Hnoland(find(H < -1)) = 0.;

%Volume
%volume(l) = dA * sum(sum(Hnoland)); %m3
%Perturbed volume
volume(l) = dA * sum(sum(eta .* mask)); %m3

%Momentum (actually it's only the integrated velocity)
iMomentumU(l) = dA * sum(sum(u_t .* Hnoland,2)); % kg * m / s
iMomentumV(l) = dA * sum(sum(v_t .* Hnoland,1)); % kg * m / s

%Kinetic energy
ke = .5 * rho0 * dA * (u_t.^2 + v_t.^2) .* Hnoland .* mask; % kg * m * m / s2 = J
iKe(l) = sum(sum(ke)); % kg * m * m / s2 = J

%Potential energy
pe = .5 * rho0 * g * dA .* eta.^2 .* mask; %kg m * m / s2 = J
iPe(l) = sum(sum(pe)); %kg m * m / s2 = J

%Eddy kinetic energy
eke = eke + rho0 * dA * K * 2 * dt * (gradux.^2 + gradvy.^2 + graduy.^2 + gradvx.^2) ...
    .* Hnoland .* mask;
iEke(l) = sum(sum(eke));

%global vorticity
iVort(l) = sum(sum(curl_w)) * dA; % m2 / s

%global enstrophy
iEnst(l) = sum(sum(enstrophy_w)) * dA; % m2 / s2

%global squared strain
iSqShear(l) = sum(sum(sqshearrate_w))* dA;
iSqStrech(l) = sum(sum(sqstrechrate_t)) * dA;
iSqStrain(l) = sum(sum(sqstrain_t)) * dA;

%global Okubo-Weiss
iOWeiss(l) = sum(sum(okuboweiss_t)) * dA; %<--- too large numerical error
%iOWeiss(l) = iSqStrain(l) - iEnst(l);

%time
vtime(l) = time;

diag = [...
        volume(l), ...
        iKe(l), ...
        iPe(l), ...
        iVort(l), ...
        iEnst(l), ...
        iSqStrain(l), ...
        iOWeiss(l) ...
        ];

function computeCurl_w
%function computeCurl_w % 1 / s

global curl_w;
global curl_t;
global mask;
global mask_v;
global mask_u,
global u;
global v;
global M;
global N;
global dy;
global dx;

dA = dx .* dy;

%Curl of velocity in the interior...
curl_w(2:M,2:N) = ((mask_v(2:M,2:N) .* v(2:M,2:N) ...
    - mask_v(1:M-1,2:N) .* v(1:M-1,2:N)) .* dy ...
    - (mask_u(2:M,2:N) .* u(2:M,2:N) ...
    - mask_u(2:M,1:N-1) .* u(2:M,1:N-1)) .* dx) ./ dA; % 1 / s
%... at the corners %Thank you for the conversation with F.L.!
%lower-left SW
curl_w(1,1) = ((mask_v(1,1) * v(1,1)) * dy ...
    - (mask_u(1,1) * u(1,1)) * dx)/dA; % 1 / s
%lower-right SE
curl_w(1,N+1) = ((mask_v(1,N+1) * v(1,N+1)) * dy ...
    - ( - mask_u(1,N) * u(1,N)) * dx)/dA; % 1 / s
%top-left NW
curl_w(M+1,1) = (( - mask_v(M,1) * v(M,1)) * dy ...
    - (mask_u(M+1,1) * u(M+1,1)) * dx)/dA; % 1 / s
%top-right NE
curl_w(M+1,N+1) = ((- mask_v(M,N+1) * v(M,N+1)) * dy ...
    - (- mask_u(M+1,N) * u(M+1,N)) * dx)/dA; % 1 / s
%... at the western boundary
curl_w(2:M,1) = ((mask_v(2:M,1) .* v(2:M,1) - mask_v(1:M-1,1) .* v(1:M-1,1)) * dy ...
    - (mask_u(2:M,1) .* u(2:M,1)) * dx)/dA; % 1 / s
%... at the eastern boundary
curl_w(2:M,N+1) = ((mask_v(2:M,N+1) .* v(2:M,N+1) - mask_v(1:M-1,N+1) .* v(1:M-1,N+1)) * dy ...
    - (- mask_u(2:M,N) .* u(2:M,N)) * dx)/dA; % 1 / s
%... at the southern boundary
curl_w(1,2:N) = ((mask_v(1,2:N) .* v(1,2:N)) * dy ...
    - (mask_u(1,2:N) .* u(1,2:N) - mask_u(1,1:N-1) .* u(1,1:N-1)) * dx)/dA; % 1 / s
%... at the northern boundary
curl_w(M+1,2:N) = ((- mask_v(M,2:N) .* v(M,2:N)) * dy ...
    - (mask_u(M+1,2:N) .* u(M+1,2:N) - mask_u(M+1,1:N-1) .* u(M+1,1:N-1)) * dx)/dA; % 1 / s
% interpolation from the {Z,W,F}-cell to the T-Cell
curl_t = interpolFtoT(curl_w, mask, M, N);

function computeShearRate_w
%function computeShearRate_w 1 / s
%

global shearrate_w;
global shearrate_t;
global mask;
global u;
global v;
global M;
global N;
global dy;
global dx;

dA = dx * dy;

%Curl of (-u,v,0)
%interior of the domain...
shearrate_w(2:M,2:N) = ((v(2:M,2:N) - v(1:M-1,2:N)) * dy ...
    + (u(2:M,2:N) - u(2:M,1:N-1)) * dx) / dA;
%bottom-left corner (SW)
shearrate_w(1,1) = ((v(1,1)) * dy ...
    + (u(1,1)) * dx) / dA;
%bottom-right corner (SE)
shearrate_w(1,N+1) = ((v(1,N+1)) * dy ...
    + (- u(1,N)) * dx) / dA;
%top-left corner (NW)
shearrate_w(M+1,1) = ((- v(M,1)) * dy ...
    + (u(M+1,1)) * dx) / dA;
%top-right corner (NE)
shearrate_w(M+1,N+1) = ((- v(M,N+1)) * dy ...
    + (- u(M+1,N)) * dx) / dA;
%western boundary
shearrate_w(2:M,1) = ((v(2:M,1) - v(1:M-1,1)) * dy ...
    + (u(2:M,1)) * dx) / dA;
%eastern boundary
shearrate_w(2:M,N+1) = ((v(2:M,N+1) - v(1:M-1,N)) * dy ...
    + (- u(2:M,N)) * dx) / dA;
%southern boundary
shearrate_w(1,2:N) = ((v(1,2:N)) * dy ...
    + (u(1,2:N) - u(1,1:N-1)) * dx) / dA;
%northern boundary
shearrate_w(M+1,2:N) = ((- v(M,2:N)) * dy ...
    + (u(M+1,2:N) - u(M+1,1:N-1)) * dx) / dA;
% interpolate from F-cell to a T-cell
shearrate_t = interpolFtoT(shearrate_w, mask, M, N);

function computePotentialVorticity_w
%function computePotentialVorticity_w % 1 / s / m
%

global potvorticity_t;
global mask;
global curl_w;
global potvorticity_w;
global coriolis;
global f;
global H;
global M;
global N;

%Adding coriolis
if coriolis
  potvorticity_w = curl_w + f; % ??? curl_w * dA  + f???
end

%potential vorticity
%within the domain...
potvorticity_w(2:M,2:N) = potvorticity_w(2:M,2:N) * 4. ./ ...
    (H(2:M,2:N) + H(1:M-1,2:N) + H(1:M-1,1:N-1) + H(2:M,1:N-1));
%bottom-left corner (SW)
potvorticity_w(1,1) = potvorticity_w(1,1) / ...
    (H(1,1));
%bottom-right corner (SE)
potvorticity_w(1,N+1) = potvorticity_w(1,N+1)/ ...
    (H(1,N));
%top-left corner (NW)
potvorticity_w(M+1,1) = potvorticity_w(M+1,1) / ...
    (H(M,1));
%top-right corner (NE)
potvorticity_w(M+1,N+1) = potvorticity_w(M+1,N+1) / ...
    (H(M,1));
%western boundary
potvorticity_w(2:M,1) = potvorticity_w(2:M,1) * 2. ./ ...
    (H(2:M,1) + H(1:M-1,1));
%eastern boundary
potvorticity_w(2:M,N+1) = potvorticity_w(2:M,N+1) * 2. ./ ...
    (H(1:M-1,N) + H(2:M,N));
%southern boundary
potvorticity_w(1,2:N) = potvorticity_w(1,2:N) * 2. ./ ...
    (H(1,2:N) + H(1,1:N-1));
%northern boundary
potvorticity_w(M+1,2:N) = potvorticity_w(M+1,2:N) * 2. ./ ...
    (H(M,2:N) + H(M,1:N-1));
%interpolate from an F-cell to a T-cell
potvorticity_t = interpolFtoT( potvorticity_w, mask, M, N);

function tgrid = interpolFtoT(fgrid, mask, M, N)
%% function tgrid = interpolFtoT(fgrid, mask, M, N)
tgrid = .25 * mask .* ( ...
              fgrid(1:M,1:N) ...
            + fgrid(2:M+1,1:N) ...
            + fgrid(2:M+1,2:N+1) ...
            + fgrid(1:M,2:N+1) ...
            );

%--------------------------------------------------------------------------
%SHEL SHallow-water numerical modEL
%Copyright (C) 2006,2009,2010,2011. Guillaume Riflet,
%Instituto Superior Técnico da Universidade Técnica de Lisboa.
%--------------------------------------------------------------------------
        