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
%     GNU General Public License for more dr.etails.
% 
%     You should have received a copy of the GNU General Public License
%     along with SHEL.  If not, see <http://www.gnu.org/licenses/>.
%--------------------------------------------------------------------------

function s = model_handles2(r)
%function s = model_handles2
%
%Call it this way:
% s = model_handles;
% y = s.model_compute_input(); % --> the parenthesis ARE REQUIRED!
% y = s.model_compute_output( handles );

%globals are replaced by a big structure
r.time              = nan;
r.timetr            = nan;
r.dt                = nan;
r.dttr              = nan;
r.ltr               = nan;
r.frame             = nan;
r.duration          = nan;
r.outputdt          = nan;
r.L                 = nan;
r.outputL           = nan;

r.frame             = nan;
r.printG            = nan;
r.film              = nan;
r.mov               = nan;
r.myfullfile        = nan;
r.mydir             = nan;
r.optprint          = nan;

r.dx                = nan;
r.dy                = nan;
r.M                 = nan;
r.N                 = nan;
r.x0                = nan;
r.y0                = nan;
r.x                 = nan;
r.y                 = nan;
r.x_u               = nan;
r.y_u               = nan;
r.x_v               = nan;
r.y_v               = nan;
r.x_w               = nan;
r.y_w               = nan;

r.loadbathymetry    = nan;
r.bathymetryfile    = nan;
r.varsinfile        = nan;

%Test-cases
r.tc_bump           = nan;
r.tc_geostrophic    = nan;
r.tc_isla           = nan;
r.tc_taylor         = nan;

%boundary condition
r.u_closed          = nan;
r.v_closed          = nan;
%r.bc_closed;
r.noslip            = nan;

%r.step bathym
r.step              = nan;

%bump
r.bump_x0           = nan;
r.bump_y0           = nan;
r.bump_sx           = nan;
r.bump_sy           = nan;
r.bump_d0           = nan;

%island
r.isla_x0           = nan;
r.isla_y0           = nan;
r.isla_R            = nan;

%r.tracer
r.tracer            = nan;
r.TrXo              = nan;
r.TrYo              = nan;
r.TrR               = nan;

%Grid and bathymetry related
r.land              = nan;
r.d0                = nan;
r.d0_r              = nan;

%T-cells
r.Tr                = nan;
r.eta0              = nan;
r.eta_old           = nan;
r.eta               = nan;
r.eta_new           = nan;
r.H_old             = nan;
r.H                 = nan;
r.H_new             = nan;
r.d                 = nan;
r.mask              = nan;
r.masknan           = nan;
r.u_t               = nan;
r.v_t               = nan;

 %Diagnostic quantities
r.ke                = nan;
r.pe                = nan;
r.eke               = nan;
r.gradux            = nan;
r.graduy            = nan;
r.gradvx            = nan;
r.gradvy            = nan;
r.curl_t            = nan;
r.potvorticity_t    = nan;
r.strechrate_t      = nan;
r.shearrate_t       = nan;
r.sqstrechrate_t    = nan;
r.sqshearrate_t     = nan;
r.enstrophy_t       = nan;
r.sqstrain_t        = nan;
r.divergence_t      = nan;
r.sqdivergence_t    = nan;
r.okuboweiss_t      = nan;

%r.u-cells
r.u_a               = nan;
r.u_old             = nan;
r.u                 = nan;
r.u_new             = nan;
r.mask_u            = nan;

%r.v-cells
r.v_a               = nan;
r.v_old             = nan;
r.v                 = nan;
r.v_new             = nan;
r.mask_v            = nan;

%Z-cells
r.curl_w            = nan;
r.potvorticity_w    = nan;
r.shearrate_w       = nan;
r.sqshearrate_w     = nan;
r.enstrophy_w       = nan;
r.okuboweiss_w      = nan;

%r.time-dependent r.properties (domain integrated)
r.volume            = nan;
r.iKe               = nan;
r.iPe               = nan;
r.vtime             = nan;
r.iVort             = nan;
r.iEnst             = nan;
r.iSqStrech         = nan;
r.iSqShear          = nan;
r.iSqStrain         = nan;
r.iOWeiss           = nan;
r.iMomentumU        = nan;
r.iMomentumV        = nan;
r.iEke              = nan;

r.K                 = nan;
r.gama              = nan;

s.model_compute_input = @ComputeModelInput;
s.model_compute_output = @ComputeModelOutput;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Here starts the model Inputs part
%%Initial conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ComputeModelInput(r)
%function ComputeModelInputs

    %GR : clears figure of images
    r.frame = 1;
    r.time = 0.;
    r.timetr = r.time + 2 * r.dt;
    ComputeTime(r);
    makecoordinates(r);
    fillfields(r);

function ComputeTime(r)
%function updateTime

    r.L = getL(r.duration,r.dt);
    r.outputL = getL(r.outputdt, r.dt);

function L_L = getL (duration_L, dt_L)
%function L_L = getL (duration_L, dt_L)

    L_L = floor(duration_L / dt_L);
    
function makecoordinates(r)
%function makecoordinates

if r.loadbathymetry
    r.varsinfile = who('-file',r.bathymetryfile);
    if ~isempty(r.varsinfile)
        load(r.bathymetryfile, 'r.d');
        siz = size(r.d);
        r.M = siz(1);
        r.N = siz(2);
    end
end

r.x0 = r.dx/2;
r.y0 = r.dy/2;

%T-cell
r.x = (1:r.M)' * (1:r.N);
r.y = (1:r.M)' * (1:r.N);
for j = 1:r.N
    r.x(:,j) = r.x0 + (r.x(:,j)/j - 1) * r.dx;
end
for i = 1:r.M
    r.y(i,:) = r.y0 + (r.y(i,:)/i - 1) * r.dy;
end

%r.u-cell
r.x_u = (1:r.M+1)' * (1:r.N);
r.y_u = (1:r.M+1)' * (1:r.N);
for j = 1:r.N
    r.x_u(:,j) = r.x0 + (r.x_u(:,j)/j - 1) * r.dx - r.dx/2;
end
for i = 1:r.M
    r.y_u(i,:) = r.y0 + (r.y_u(i,:)/i - 1) * r.dy;
end

%r.v-cell
r.x_v = (1:r.M)' * (1:r.N+1);
r.y_v = (1:r.M)' * (1:r.N+1);
for j = 1:r.N
    r.x_v(:,j) = r.x0 + (r.x_v(:,j)/j - 1) * r.dx;
end
for i = 1:r.M
    r.y_v(i,:) = r.y0 + (r.y_v(i,:)/i - 1) * r.dy - r.dy/2;
end

%W-cell
r.x_w = (1:r.M)' * (1:r.N);
r.y_w = (1:r.M)' * (1:r.N);
for j = 1:r.N
    r.x_w(:,j) = r.x0 + (r.x_w(:,j)/j - 1) * r.dx - r.dx/2;
end
for i = 1:r.M
    r.y_w(i,:) = r.y0 + (r.y_w(i,:)/i - 1) * r.dy - r.dy/2;
end

function fillfields(r)
%function fillfields

%         Arakawa C staggered grid (Arakawa 1966)
%
%    --j y
%   |           --- r.v ---           r.v ------- r.v         r.u -- T -- r.u
%   i          |         |          |         |         |         |
%   x          r.u    T    r.u          T    r.u    T         |    r.v    |
%              |         |          |         |         |         |
%               --- r.v ---           r.v ------- r.v         r.u -- T -- r.u
%                 T-cell               r.u-cell              r.v-cell
%
%       T(r.eta, r.H, r.Tr, r.d, f, r.mask,   r.u(r.u, r.mask_u         r.v(r.v, r.mask_v
%           x, y, wind and bottom      r.x_u, r.y_u)           r.x_v, r.y_v)
%            stress)
%    --j y
%   |          T--- r.u ---T 
%   i          |         | 
%   x          r.v    W    r.v 
%              |         | 
%              T--- r.u ---T 
%                 W-cell   
%
%       W(zr.eta - curl of r.v, vorticity)


%Depth field (bathymetry)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r.land = -99;

%%Constant depth or read from bathymetry file
if r.tc_taylor %Submarine mount equals 30% of total mean depth
    r.d = r.d0 - cylindermap(r,zeros(size(r.d)), 0.3 * r.d0, r.bump_sx, r.bump_sy, ...
                            r.M * r.dx * .5 - r.bump_x0, ...
                            r.N * r.dy * .5 - r.bump_y0);
elseif r.loadbathymetry && sum(strcmp(r.varsinfile,'r.d'))
    load(r.bathymetryfile, 'r.d');
else
    r.d = r.d0 * ones(r.M,r.N);
end

%%changing depth
if r.step
    r.d = maker.step(r.d,r.d0_r.step);
end

%%Closed boundaries: boundaries are filled with r.land
if r.v_closed
    for j = 1:r.N-1:r.N;
    for i = 1:r.M
        r.d(i,j) = r.land;
    end
    end
end

if r.u_closed    
    for i = 1:r.M-1:r.M;
    for j = 1:r.N
        r.d(i,j) = r.land;
    end
    end
end

%%Inner islands (random)
%%%arguments: depth, io, jo, r
if r.tc_isla
    r.d = makeisland(r,r.d, r.M * r.dx * .5 - r.isla_x0, r.N * r.dy * .5 - r.isla_y0, r.isla_R);
end

%WATERMASK matrices%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%T-cell
%%all-water
r.mask = ones(r.M,r.N);

%%eventual islands and closed domain
for j=1:r.N
for i=1:r.M
    if r.d(i,j) < 0.
        r.mask(i,j) = 0;
    end
end
end

%r.u-Cell and r.v-cell
r.mask_u = ones(r.M+1,r.N);
r.mask_v = ones(r.M,r.N+1);
for j = 1:r.N
for i = 1:r.M
    if r.mask( i, j) == 0
        r.mask_u(  i,  j) = r.mask(i,j);
        r.mask_u(i+1,  j) = r.mask(i,j);
        r.mask_v(  i,  j) = r.mask(i,j);
        r.mask_v(  i,j+1) = r.mask(i,j);
    end
end
end

if r.noslip
    %r.u-cells
    for j = 2:r.N-1
    for i = 2:r.M
        r.mask_u(i,j) = r.mask(i  ,j+1) * r.mask(i  ,j-1) * r.mask(i-1,j+1) * r.mask(i-1,j-1);
    end
    end
    %r.v-Cells
    for j = 2:r.N
    for i = 2:r.M-1
        r.mask_v(i,j) = r.mask(i+1,  j) * r.mask(i+1,j-1) * r.mask(i-1,  j) * r.mask(i-1,j-1);
    end
    end
end

%Water level fields%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (r.tc_bump || r.tc_geostrophic)
    r.eta_old = bumpmap(r,r.bump_d0, r.bump_sx, r.bump_sy, r.M * r.dx * .5 - r.bump_x0, ...
                                                r.N * r.dy * .5 - r.bump_y0);
else
    if r.loadbathymetry && sum(strcmp(r.varsinfile,'r.eta'))
        load(r.bathymetryfile, 'r.eta');
        r.eta_old = r.eta;
    else
        r.eta_old = r.eta0 * ones(r.M,r.N);
    end
end
r.eta = r.eta_old;
r.eta_new = r.eta;

%Bottom-to-surface depth (always positive below the water surface)
%r.H = r.eta + r.d.
r.H_old = r.eta_old + r.d;
r.H =r.H_old;
r.H_new = r.H_old;

if r.tracer
    %Unitary concentration...
    r.Tr = zeros(r.M,r.N);
    r.Tr = cylindermap(r,r.Tr, 1., r.TrR, r.TrR, r.M * r.dx * .5 - r.TrXo, ...
                                    r.N * r.dy * .5 - r.TrYo);
end

%r.u velocity fields%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if r.tc_geostrophic
    r.u_old = ugeo(r.eta);    
elseif r.tc_taylor % 3 cm/s r.u-flow
    r.u_old = zeros(size(r.mask_u));
    %%uncomment one of the following lines
    r.u_old = 0.03 * ones(size(r.mask_u));
    %r.u_old(1:r.M:r.M+1,:) = 0.03 * ones(size(r.mask_u(1:r.M:r.M+1,:)));
else
    if r.loadbathymetry && sum(strcmp(r.varsinfile,'r.u'))
        load(r.bathymetryfile, 'r.u');
        r.u_old = r.u;
    else
        r.u_old = zeros(r.M+1,r.N);
    end
end
r.u = r.u_old;
r.u_new = r.u;
r.u_a = r.u;

%r.v velocity fields%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if r.tc_geostrophic
    r.v_old = vgeo(r.eta);    
else
    if r.loadbathymetry && sum(strcmp(r.varsinfile,'r.v'))
        load(r.bathymetryfile, 'r.v');
        r.v_old = r.v;
    else
        r.v_old = zeros(r.M,r.N+1);
    end
end
r.v = r.v_old;
r.v_new = r.v;
r.v_a = r.v;

%call this function when you need 
%to manually initialize the level
%in the Taylor column flow.
if r.tc_taylor
    makeGeostLevelFromU(r)
end

%Z field%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r.curl_w = zeros(r.M+1,r.N+1);
r.potvorticity_w = zeros(r.M+1,r.N+1);
r.enstrophy_w = zeros(r.M+1,r.N+1);
r.shearrate_w = zeros(r.M+1,r.N+1);
r.sqshearrate_w = zeros(r.M+1,r.N+1);
r.okuboweiss_w = zeros(r.M+1,r.N+1);

%Visualization matrices%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r.masknan = ones(r.M,r.N);
for j=1:r.N
for i=1:r.M
        if r.mask(i,j) == 0
            r.masknan(i,j) = nan;
        else
            r.masknan(i,j) = 1;
        end
end
end

r.u_t = zeros(r.M,r.N);
r.v_t = zeros(r.M,r.N);

r.ke = zeros(r.M,r.N);
r.pe = zeros(r.M,r.N);

r.curl_t = zeros(r.M,r.N);
r.potvorticity_t = zeros(r.M,r.N);
r.strechrate_t = zeros(r.M,r.N);
r.shearrate_t = zeros(r.M,r.N);
r.enstrophy_t = zeros(r.M,r.N);
r.sqstrain_t = zeros(r.M,r.N);
r.okuboweiss_t = zeros(r.M,r.N);
r.sqshearrate_t = zeros(r.M,r.N);
r.sqstrechrate_t = zeros(r.M,r.N);
r.divergence_t = zeros(r.M,r.N);
r.sqdivergence_t = zeros(r.M,r.N);

%Eddy kinetic energy
r.gradux=zeros(r.M,r.N);
r.graduy=zeros(r.M,r.N);
r.gradvx=zeros(r.M,r.N);
r.gradvy=zeros(r.M,r.N);
r.eke=zeros(r.M,r.N);

%One-dimensional diagnostic variables
%Energy initialization
r.vtime = nan(1,L);
r.iKe = nan(1,L);
r.iPe = nan(1,L);
r.iEke = nan(1,L);
r.volume = nan(1,L);
r.iVort = nan(1,L);
r.iEnst = nan(1,L);
r.iSqShear = nan(1,L);
r.iSqStrech = nan(1,L);
r.iSqStrain = nan(1,L);
r.iOWeiss = nan(1,L);
r.iMomentumU = nan(1,L);
r.iMomentumV = nan(1,L);

ComputeDiagnostics(1);

%cylinder field
function map_L = cylindermap(r,map_L, l, sx, sy, x_L, y_L)
       
    ind = find ( (x - x_L).^2 / sx^2 + (y - y_L).^2 / sy^2 < 1. );
    map_L(ind) = l;
    
%bump test-case field
function map_L = bumpmap(r,l, sx, sy, x_L, y_L)
       
    map_L = zeros(r.M,r.N);   
    
    for i = 1:r.M
    for j = 1:r.N
        map_L(i,j) = l * exp( - ( ( x(i,j) - x_L )^2 / sx^2 + ( y(i,j) - y_L )^2 / sy^2 ) );
    end
    end
        
function depth = makeisland(r,depth, x_L, y_L, r)
%function depth = makeisland(depth, x_L, y_L, r)

    for i = 1:r.M
    for j = 1:r.N
        if r^2 > ( x(i,j) - x_L )^2 + ( y(i,j) - y_L )^2
            depth(i,j) = r.land;
        end
    end
    end

function ugeo_L = ugeo(r, r.eta_L)
%function ugeo_L = ugeo(r.eta_L)

    j_L = floor(r.N/2)+1;
    sig_L = 1.;
    ug_L = zeros(r.M,r.N);
    
    for i=1:r.M
    for j=1:r.N
        ug_L(i,j) = g / f * ( j - j_L ) / sig_L^2 * r.eta_L(i,j);
    end
    end
    ugeo_L = InterpolTU(ug_L);

function vgeo_L = vgeo(r,r.eta_L)
%function vgeo_L = vgeo(r.eta_L)

    i_L = floor(r.M/2)+1;
    sig_L = 1.;
    vg_L = zeros(r.M,r.N);
    
    for i=1:r.M
    for j=1:r.N
        vg_L(i,j) = - g / f * ( i - i_L ) / sig_L^2 * r.eta_L(i,j);
    end
    end
    vgeo_L = InterpolTV(vg_L);

function bath = makestep(bath,depth)
%function bath = makestep(bath,depth)
%Creates a r.stepped bathymetry
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

function makeGeostLevelFromU(r)

    for j=2:r.N
    for i=1:r.M
        r.eta(i,j) = r.eta(i,j-1) * r.mask(i,j-1) - r.dy * f/g * (...
            r.mask_u(i,j) * r.u(i,j) + r.mask_u(i+1,j) * r.u(i+1,j) ...
            + r.mask_u(i,j-1) * r.u(i,j-1) + r.mask_u(i+1,j-1) * r.u(i+1,j-1) ...
            ) / ( ...
            r.mask_u(i,j) + r.mask_u(i+1,j) + r.mask_u(i,j-1) + r.mask_u(i+1,j-1) ...
            );
        r.eta(i,j) = r.mask(i,j) * r.eta(i,j);
    end
    end

    r.eta_old = r.eta;
    r.eta_new = r.eta;

    r.H_old = r.eta_old + r.d;
    r.H =r.H_old;
    r.H_new = r.H_old;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Here starts the model Outputs part
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ComputeModelOutput( r, handles )
%function ComputeModelOutputs( handles ) [old ComputeModel]

%r.time parameters and output parameters

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

% Main loop in r.time
r.ltr = 0;
for l = 1 : L           
    
    %Takes care of the numerical scheme
    %for momentum and for the r.tracer
    %This line can easily be replaced
    %with a switch, choosing one kernel
    %among many (mex file kernel, other
    %kernels ...).
    kernel_SHEL(r,l);
    
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

        %Display r.time
        strin = sprintf('r.time = %0.3g s', r.time);
        set(handles.timetext,'String',strin);
        pause(0.01); %s pause (required if I want to update the clock).
        
    end
    ol=ol+1;
    
%End of main loop in r.time    
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
disp(sprintf('Elapsed r.time: %0.5f s\nTime per iteration: %0.5f s\r.N', cput, cput/L));
if r.tracer
    disp(sprintf('Made %0.3d r.tracer iterations on a total of %0.3d iterations', r.ltr, l));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Here starts the Kernel part
%%the numerical scheme
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function kernel_SHEL(r,l)
%function kernel_SHEL (old ComputeModel)

%Explicit leapfrog scheme
%combined with Asselin-Roberts filter
%(the filter supresses the computational mode
%and prevents decoupling between odd/even timer.steps)
%
%         Arakawa C staggered grid (Arakawa 1966)
%
%    --j y
%   |           --- r.v ---           r.v ------- r.v         r.u -- T -- r.u
%   i          |         |          |         |         |         |
%   x          r.u    T    r.u          T    r.u    T         |    r.v    |
%              |         |          |         |         |         |
%               --- r.v ---           r.v ------- r.v         r.u -- T -- r.u
%                 T-cell               r.u-cell              r.v-cell
%                  r.M x r.N               r.M+1 x r.N             r.M x r.N+1
%       T(r.eta, r.H, r.Tr, r.d, f, r.mask,   r.u(r.u, r.mask_u         r.v(r.v, r.mask_v
%           x, y, wind and bottom      r.x_u, r.y_u)           r.x_v, r.y_v)
%            stress, curl(r.v,r.u,0))
%
%    --j y
%   |           --- r.u ---  
%   i          |         | 
%   x          r.v    W    r.v 
%              |         | 
%               --- r.u ---  
%                 T-cell
%                r.M+1 x r.N+1
%       W( curl(r.u,r.v,0) )
%        (deformation rate = curl(-r.u,r.v,0)
%%%
   
    %Update the r.time
    r.time = r.time + r.dt;
    
    %Compute the r.time scheme of the momentum+continuity model
    ComputeLeapfrog;

    %Compute all the tracers
    if r.tracer
        if r.timetr <= r.time
            [r.Tr, r.dttr] = ComputeTracer_FT(r.Tr, r.K);
            r.timetr = r.time + r.dttr;
            r.ltr = r.ltr + 1;
        end
    end
    
    %Numerics to compute the diagnostic variables
    ComputeDiagnostics(r, l);

function ComputeLeapfrog(r)
%Shallow-water equations solver
%
 
%T-cells
r.dtmask;
 
%r.u-cells
r.dtmask_u;
 
%r.v-cells
r.dtmask_v;

%% Solve the spatial schemes in the interior of the domain
    % For the water elevation
    RHSr.eta = ComputeContinuity(r.H, r.eta, r.u, r.v, r.mask);
    r.eta_new = r.eta_old + 2 * r.dt * RHSr.eta .* r.mask; 
    r.H_new = r.eta_new + r.d;
    ind = find(r.H_new < 0);
    r.eta_new(ind) = -r.d;
    dtmask = (r.eta_new - r.eta_old) * .5 ./ RHSr.eta .* r.mask;
    %Must build the dtmask_u and dtmask_v in order to preserve momentum
    %too.
    dtmask_u = buildmask_U(dtmask,r.mask_u);
    dtmask_v = buildmask_V(dtmask',r.mask_v')';
    
    RHSr.eta = ComputeContinuity(r.H, r.eta, r.u, r.v, dtmask);
    
    % For the r.u and r.v components of velocity
    % % The old 10%-faster way (but a lot more error-prone)
    %     [RHSu, RHSv] = ComputeUV2D_old;
    % Switch r.v <-> r.u; r.v_old <-> r.u_old; r.dy <-> r.dx and -1. <-> 1. (coriolis)
    % Transpose the matrices.
    RHSu = ComputeSpaceU_CS( r.H, r.H_old, r.eta, ...
                        r.u, r.u_old, r.v, r.v_old, ...
                        dtmask_u, r.dx, r.dy, 1.);
    RHSv = ComputeSpaceU_CS( r.H', r.H_old', r.eta', ...
                        r.v', r.v_old', r.u', r.u_old', ...
                        dtmask_v', r.dy, r.dx, -1.)';
 
%% LEAPFROG r.time-SCHEME
    % Water level r.time scheme
    r.eta_new = r.eta_old + 2 * RHSr.eta .* r.mask; 
    r.H_new = r.eta_new + r.d;
    %r.u and r.v r.time schemes inside the domain. 
    %Leapfrog means that 2xdt is passed as teh r.time r.step.
    r.v_new = ComputeTimeU_FT(r.v_old', r.v_new', RHSv', r.H_old', r.H_new', r.mask_v', 2.)';
    r.u_new = ComputeTimeU_FT(r.u_old, r.u_new, RHSu, r.H_old, r.H_new, r.mask_u, 2.);

%% r.time scheme at the northern/southern and eastern/southern OB (2 x r.d
%% <- LeapFrog)

    %Southern and northern boundaries (along r.v -> transpose inputs)
    %(transposed inputs -> transposed outputs).
    [r.eta_new, r.H_new, r.v_new, r.u_new] = ComputeOB_U( ...
            r.eta_new', r.eta_old', r.H_new', r.H_old', r.d',...
            r.v_new', r.mask_v', ...
            r.u_new', r.u_old', r.mask_u', ...
            2*r.dt, r.dy ...
     );
    %Eastern and western boundaries (along r.u)
    %The previous outputs were transposed, so we transpose
    %them back as inputs.
    [r.eta_new, r.H_new, r.u_new, r.v_new] = ComputeOB_U( ...
            r.eta_new', r.eta_old, r.H_new', r.H_old, r.d,...
            r.u_new', r.mask_u, ...
            r.v_new', r.v_old, r.mask_v, ...
            2*r.dt, r.dx ...
     );
     
%% Asselin-Robert filter 
    r.eta = r.eta + r.gama * ( r.eta_old - 2 * r.eta + r.eta_new ); 
    r.u = r.u + r.gama * ( r.u_old - 2 * r.u + r.u_new ); 
    r.v = r.v + r.gama * ( r.v_old - 2 * r.v + r.v_new ); 

%% Update matrices
    r.eta_old = r.eta;
    r.eta = r.eta_new;
     
    r.H_old = r.H;
    r.H = r.H_new;
 
    r.u_old = r.u;
    r.u = r.u_new;
     
    r.v_old = r.v;
    r.v = r.v_new;

function RHSr.eta = ComputeContinuity(r.H, r.eta, r.u, r.v, r.mask)
%function ComputeContinuity
%Shallow-water equations solver

% r.H; 
% r.eta; 
% r.u; 
% r.v; 
% r.mask; 
% rHSr.eta;
r.dx; 
r.dy; 

[r.M,r.N] = size(r.eta);
RHSr.eta = zeros(r.M,r.N); 

%Compute along r.u then along r.v
RHSr.eta(2:r.M-1,2:r.N-1) = - ( ...
    ComputeContinuityU(r.H,r.u,r.mask,r.dx) + ComputeContinuityU(r.H',r.v',r.mask',r.dy)' ...
);

function RHS = ComputeContinuityU(r.H,r.u,r.mask,r.dx)
[r.M,r.N] = size(r.H);
RHS = ( ... 
            r.mask(3:r.M,2:r.N-1) .* .5 ...
            .* ( r.H(3:r.M,2:r.N-1) + r.H(2:r.M-1,2:r.N-1) ) .* r.u(3:r.M,2:r.N-1) ... %i+1, i, i+1 
            - r.mask(1:r.M-2,2:r.N-1) .* .5 ...
            .* ( r.H(2:r.M-1,2:r.N-1) + r.H(1:r.M-2,2:r.N-1) ) .* r.u(2:r.M-1,2:r.N-1) ... %i, i-1, i 
        )/r.dx;

function RHSu = ComputeSpaceU_CS( ...
                    r.H, r.H_old, r.eta, ...
                    r.u, r.u_old, r.v, r.v_old, ...
                    r.mask_u, r.dx, r.dy, ...
                    signcoriolis ...
                )
%% function newu = ComputeSpaceU_CS
%Shallow-water equations solver for the r.u component
%of the flow Centred in Space (CS).
%
% Centred differences (CD) or Centred in Space (CS).
%
%NOTE: Any viscous and frictional term
%must be evaluated in tindex(i-1,r.M)e l-1.
%
% switch r.v <-> r.u; j <-> i and r.dy <-> r.dx relative to ComputeU2D
% switch r.v_old <-> r.u_old
% take special care for the coriolis force term

%% r.variables
r.g; 
r.f; 
r.Nu; 
r.wind; 
r.bottom; 
r.pressure; 
r.coriolis; 
r.wall;

r.uwind; 
r.vwind; 
rho0;
rho_air;
r.lb; 
r.karman; 

wall = false; 

%% r.u spatial scheme
[r.M, r.N] = size(r.eta);
RHSu = zeros(size(r.u));
r.H_u = zeros(size(r.u));
r.H_old_u = zeros(size(r.u));
v_old_u = zeros(size(r.u));
v_u = zeros(size(r.u));
r.H_u(2:r.M,2:r.N-1) = twoaverage_u(r.H,r.M,r.N);
r.H_old_u(2:r.M,2:r.N-1) = twoaverage_u(r.H_old,r.M,r.N);
v_old_u(2:r.M,2:r.N-1) = fouraverage_u(r.v_old(:,2:r.N),r.M, r.N-2); 
v_u(2:r.M,2:r.N-1) = fouraverage_u(r.v(:,2:r.N),r.M, r.N-2); 

RHSu(2:r.M,2:r.N-1) = r.mask_u(2:r.M,2:r.N-1) .* ( ... 
    ...%advection along r.u ... 
        - 0.2500 / r.dx * ( ... %i-1, i 
        r.mask_u(3:r.M+1,2:r.N-1) .* ( r.u(3:r.M+1,2:r.N-1) + r.u(2:r.M,2:r.N-1) ).^2 .* r.H(2:r.M,2:r.N-1) ... %i+1, i, i+1, i, i 
      - r.mask_u(1:r.M-1,2:r.N-1) .* ( r.u(2:r.M,2:r.N-1) + r.u(1:r.M-1,2:r.N-1) ).^2 .* r.H(1:r.M-1,2:r.N-1) ... %i, i-1, i, i-1, i-1 
     )  ... 
    ...%advection along r.v 
     - 0.0625 / r.dy * ( ... %i-1, i 
        ( ... 
            r.mask_u(2:r.M,3:r.N) .* ( r.H(2:r.M,3:r.N) + r.H(2:r.M,2:r.N-1) + r.H(1:r.M-1,3:r.N) + r.H(1:r.M-1,2:r.N-1) ) ... % j+1 ; j ; i-1,j+1 ; i-1 
          .* ( r.v(2:r.M,3:r.N) + r.v(1:r.M-1,3:r.N) ) .* ( r.u(2:r.M,3:r.N) + r.u(2:r.M,2:r.N-1) ) ...% j+1 ; i-1,j+1 ; j+1 ; i 
          - r.mask_u(2:r.M,1:r.N-2) .* ( r.H(2:r.M,2:r.N-1) + r.H(2:r.M,1:r.N-2) + r.H(1:r.M-1,2:r.N-1) + r.H(1:r.M-1,1:r.N-2) ) ... % i ; j-1 ; i-1 ; i-1,j-1 
          .* ( r.v(2:r.M,2:r.N-1) + r.v(1:r.M-1,2:r.N-1) ) .* ( r.u(2:r.M,2:r.N-1) + r.u(2:r.M,1:r.N-2) ) ...% i ; i-1,j-1 ; i ; j-1 
        ) ... 
     ) ... 
    ...%Coriolis Force (r.u -> +; r.v -> -)
    + coriolis * ( ... 
        signcoriolis * f * v_u(2:r.M,2:r.N-1) .* r.H_u(2:r.M,2:r.N-1) ... % i ; i-1 . 
    )... 
    ...%Pressure gradient Force 
    + pressure * ( ... 
       - g * r.H_u(2:r.M,2:r.N-1) .* ( r.eta(2:r.M,2:r.N-1) - r.eta(1:r.M-1,2:r.N-1) ) / r.dx ...  
                               ... % i ; i-1 ;                 i ; i-1 
    )... 
    ...%Viscosity along r.u 
      + ( ... % i ; i-1 
          + Nu * r.mask_u(3:r.M+1,2:r.N-1) .* r.H_old(2:r.M,2:r.N-1) .* ( r.u_old(3:r.M+1,2:r.N-1) - r.u_old(2:r.M,2:r.N-1) ) / r.dx ... % i ; i+1 ; i ; i 
          - Nu * r.mask_u(1:r.M-1,2:r.N-1) .* r.H_old(1:r.M-1,2:r.N-1) .* ( r.u_old(2:r.M,2:r.N-1) - r.u_old(1:r.M-1,2:r.N-1) ) / r.dx ... % i-1 ; i ; i-1; i-1 
      ) / r.dx ... % i ; i-1 
    ...%Viscosity along r.v 
      + ( ... % i ; i-1 ; i ; i-1 . 
          + Nu / r.dy * r.mask_u(2:r.M,3:r.N) .* r.H_old(2:r.M,2:r.N-1) .* ( r.u_old(2:r.M,3:r.N) - r.u_old(2:r.M,2:r.N-1) ) ... % (j;j+1;i-1,j+1;i-1);j+1 ; j 
          - Nu / r.dy * r.mask_u(2:r.M,1:r.N-2) .* r.H_old(2:r.M,1:r.N-2) .* ( r.u_old(2:r.M,2:r.N-1) - r.u_old(2:r.M,1:r.N-2) ) ... % (j-1;j;i-1;i-1,j-1);j ; j-1 
      ) / r.dy ...
);
 
%Add bottom stress. Look out, so that the bottom stress doesn't
%reverts the direction of momentum ... This wouldn't be very realistic.
if bottom || wall 
    stress = zeros(size(RHSu(2:r.M,2:r.N-1))); 
    if bottom 
        stress = bottomstress( r.u_old(2:r.M,2:r.N-1), v_old_u(2:r.M,2:r.N-1), ...
                                r.H_old_u(2:r.M,2:r.N-1), ...
                                lb, karman ); 
    end 
    if wall 
        stress = stress + bottomstress(r.u_old(2:r.M,2:r.N-1), 0, r.dy, lb, karman) ...
            .* ( ... 
                    (1 - r.mask(1:r.M-1,3:r.N)) .* r.H_old(1:r.M-1,3:r.N) + ... 
                    (1 - r.mask(2:r.M,3:r.N)) .* r.H_old(2:r.M,3:r.N) + ... 
                    (1 - r.mask(1:r.M-1,1:r.N-2)) .* r.H_old(1:r.M-1,1:r.N-2) + ... 
                    (1 - r.mask(2:r.M,1:r.N-2)) .* r.H_old(2:r.M,1:r.N-2) ... 
            ) * 0.5 / r.dy; 
    end 
    ind = find( (RHSu(2:r.M,2:r.N-1) - r.mask_u(2:r.M,2:r.N-1) .* stress) ./ RHSu(2:r.M,2:r.N-1) > 0 );
    RHSu(ind) = RHSu(ind) - r.mask_u(ind) .* stress(ind); 
    ind = find( (RHSu(2:r.M,2:r.N-1) - r.mask_u(2:r.M,2:r.N-1) .* stress) ./ RHSu(2:r.M,2:r.N-1) <= 0 );
    RHSu(ind) = 0.; 
end
    
%Add wind stress.
if wind 
    RHSu(2:r.M,2:r.N-1) = RHSu(2:r.M,2:r.N-1) ...
            + r.mask_u(2:r.M,2:r.N-1) .* windstress(...
                uwind * ones(size(r.u(2:r.M,2:r.N-1))), ...
                vwind * ones(size(r.u(2:r.M,2:r.N-1))), ...
                rho0, rho_air ...
            );
end
    
function r.u_new = ComputeTimeU_FT(r.u_old, r.u_new, RHSu, r.H_old, r.H_new, r.mask_u, r.dt)
%function r.u_new = ComputeTimeU_FT(r.u_old, r.u_new, RHSu, r.H_old, r.H_new, r.mask_u, r.dt)        
%
%Numerical solver with the Forward in r.time (FT) Euler r.time r.stepping method.
%
    [r.M,r.N] = size(r.H_old);
    
    r.u_new(2:r.M,2:r.N-1) = r.mask_u(2:r.M,2:r.N-1) .* ( ... 
                r.u_old(2:r.M,2:r.N-1) .* ( r.H_old(2:r.M,2:r.N-1) + r.H_old(1:r.M-1,2:r.N-1) ) ... 
                + 2 * r.dt * RHSu(2:r.M,2:r.N-1) ... 
              ) ... 
              ./ ( r.H_new(2:r.M,2:r.N-1) + r.H_new(1:r.M-1,2:r.N-1) );
          
function [r.eta_new, r.H_new, r.u_new, r.v_new] = ComputeOB_U( ...
    r.eta_new, r.eta_old, r.H_new, r.H_old, r.d,...
    r.u_new, r.mask_u, ...
    r.v_new, r.v_old, r.mask_v, ...
    r.dt, r.dx ...
)
%function computeOB
%Solve the north, south, east and west boundary condition
%for the r.eta, r.u and r.v variables.

r.g; 
  
%BC conditions 
r.flather;
r.neumann;
radiatr.etan; 
radiatelevel;

flather = 0;
neumann = 0;
radiatr.etan = 0; 
radiatelevel = 1;

[r.M,r.N] = size(r.d);

%% Eastern (1,:) and Western (r.M,:) boundaries    
if radiatelevel 
    %Waterlevel: GW radiation 
    r.eta_new(1:r.M-1:r.M,:) = r.mask_u(1:r.M:r.M+1,:) .* ( ...
        r.eta_old(1:r.M-1:r.M,:)- r.dt / r.dx .* sqrt( g * r.H_old(1:r.M-1:r.M,:) ) ... 
        .* ( r.eta_old(1:r.M-1:r.M,:) - r.eta_old(2:r.M-3:r.M-1,:) ) ...
    ); 
    %no r.mask below, cuz r.H must yield -99 where r.land is.W
    r.H_new(1:r.M-1:r.M,:) = r.eta_new(1:r.M-1:r.M,:) + r.d(1:r.M-1:r.M,:);
end 
if flather 
    %Normal velocity: Flather 
    r.u_new(1:r.M:r.M+1,:) = [-1 * ones(1,r.N); ones(1,r.N)] .*  r.mask_u(1:r.M:r.M+1,:) ...
        .* sqrt( g ./ r.H_new(1:r.M-1:r.M,:) ) .* r.eta_new(1:r.M-1:r.M,:); 
elseif neumann
    %Normal velocity: Null-gradient 
    r.u_new(1:r.M:r.M+1,1:r.N) = r.mask_u(1:r.M:r.M+1,1:r.N) .* r.u_new(2:r.M-2:r.M,1:r.N); 
end 
if radiatr.etan 
    %Tangential velocity: GW radiation 
    r.v_new(1:r.M-1:r.M,2:r.N) = r.mask_v(1:r.M-1:r.M,2:r.N) .* ... 
        ( ... 
            r.v_old(1:r.M-1:r.M,2:r.N) .* ( r.H_old(1:r.M-1:r.M,2:r.N) + r.H_old(1:r.M-1:r.M,1:r.N-1) ) ... 
            - r.dt / r.dx .* sqrt( g * .5 * ( r.H_old(1:r.M-1:r.M,2:r.N) + r.H_old(1:r.M-1:r.M,1:r.N-1) ) ) ... 
            .* ( r.v_old(1:r.M-1:r.M,2:r.N) - r.v_old( 2:r.M-3:r.M-1,2:r.N) ) ... 
        ) ... 
        ./ ( r.H_new(1:r.M-1:r.M,2:r.N) + r.H_new(1:r.M-1:r.M,1:r.N-1) ); 
end


function [r.Tr, r.dttr] = ComputeTracer_FT(r.Tr, K_L)
%Advection-diffusion r.tracer equation solver
%Explicit forward in r.time and upwind
%
r.dx;
r.dy;
r.dt;
 
%T-cells
r.H;
r.H_new;
r.mask;
 
%r.u-cells
r.u_a;
r.mask_u;

%r.v-cells
r.v_a;
r.mask_v;

%% Solve the r.tracer spatial scheme in the whole of the domain
    RHSt = r.mask .* ( ...
        ComputeSpaceT_UP_U(r.H,r.u_a,r.mask_u, r.Tr, K_L, r.dx) ... %Along r.u
       + ComputeSpaceT_UP_U(r.H',r.v_a',r.mask_v',r.Tr',K_L,r.dy)' ... %then along r.v
    );

%% FIND THE OPTIMAL r.tracer r.dt AS A MULTIPLE OF MOMENTUM r.dt
% within the reasonable limit of maximum 100 * r.dt.
    minRHSt  = min(min(RHSt));
    ind = find(RHSt == minRHSt);
    r.dttr = min( floor( abs(- r.H(ind) .* r.Tr(ind) ./ minRHSt) / r.dt ), 1000) ...
        * r.dt;

%% r.tracer FORWARD-IN-r.time SCHEME
    Tr_new = r.mask .* (r.H .* r.Tr + r.dttr * RHSt) ./ r.H_new;

%% Update matrices
    r.Tr = Tr_new;

function RHSt = ComputeSpaceT_UP_U(r.H, r.u, r.mask_u, r.Tr, K_L, dx_L)
%function RHSt = ComputeSpaceT_UP_U(r.H, r.u, r.mask_u, r.Tr, r.K, r.dx)
%
%UPWIND
%
    r.nullgrad;
    nullgrad = false;
    
    [r.M,r.N]=size(r.H);
    Hm = zeros(r.M+1,r.N);
    TrUm = zeros(r.M+1,r.N);
    TrUp = zeros(r.M+1,r.N);
    TrDif = zeros(r.M+1,r.N);

    %Computing quantities for the interior of the domain
    TrUm(1:r.M,:) = .5 * (abs(r.u(1:r.M,:)) - r.u(1:r.M,:)) .* r.Tr;
    TrUp(2:r.M+1,:) = .5 * (abs(r.u(2:r.M+1,:)) + r.u(2:r.M+1,:)) .* r.Tr;
    TrDif(2:r.M,:) = r.Tr(1:r.M-1,:) - r.Tr(2:r.M,:);
    Hm(2:r.M,:) = .5 * ( r.H(1:r.M-1,:) + r.H(2:r.M,:) );
    
    %Computing the quantities for the boundaries
    if nullgrad %Considering null-gradient:
        TrUm(r.M+1,:) = .5 * (abs(r.u(r.M+1,:)) - r.u(r.M+1,:)) .* r.Tr(r.M,:);
        TrUp(1,:) = .5 * (abs(r.u(1,:)) + r.u(1,:)) .* r.Tr(1,:);
        TrDif(1:r.M:r.M+1,:) = 0.;
    else %else null value is considered
        TrUm(r.M+1,:) = .5 * (abs(r.u(r.M+1,:)) - r.u(r.M+1,:)) * 0.;
        TrUp(1,:) = .5 * (abs(r.u(1,:)) + r.u(1,:)) * 0.;
        TrDif(1,:) = - r.Tr(1,:);
        TrDif(r.M+1,:) = r.Tr(r.M,:);
    end
    
    %This algo considers that the bathymetry has always a null gradient
    %at the boundary:
    Hm(1:r.M:r.M+1,:) = r.H(1:r.M-1:r.M,:);

    %Computing the fluxes at each face of the r.M+1 faces
    FFluxU = ComputeFaceFluxT_UP(Hm, TrUp, TrUm, TrDif, dx_L, K_L);
    
    %The T-cell increment is the flux balance.
    RHSt = r.mask_u(1:r.M,:) .* FFluxU(1:r.M,:) ...
         - r.mask_u(2:r.M+1,:) .* FFluxU(2:r.M+1,:);

function RHSt = ComputeFaceFluxT_UP(Hm, TrUp, TrUm, TrDif, dx_L, K_L)
%function RHSt = ComputeFaceFluxT_UP(Hm, TrUp, TrUm, TrDif, dx_L, K_L)
%Pure upwind scheme.
    RHSt = Hm .* ( ...
            TrUp - TrUm ...
            + K_L / dx_L * TrDif ...
    ) / dx_L;
        
function av = fouraverage_u(r.v, r.M, r.N) 
%% function av = fouraverage_u(r.v, r.mask, r.M, r.N)

    av =  .25 * ( ... 
            r.v(2:r.M,1:r.N) + ... 
            r.v(2:r.M,2:r.N+1) + ... 
            r.v(1:r.M-1,1:r.N) + ... 
            r.v(1:r.M-1,2:r.N+1) ... 
           ); 

function av = twoaverage_u( r.H, r.M, r.N)
    av = .5 * ( ...
            r.H(2:r.M,2:r.N-1) + r.H(1:r.M-1,2:r.N-1) ...
         );

function tb_L = bottomstress( u_L, v_L, r.H_L, lb, karman) 
% Bottom stress        
%function tb_L = bottomstress( u_L, v_L, r.H_L ) 
% 
%r.H_L is the total depth height 
     
    cd_L = karman / log( ( 0.5 .* r.H_L + lb ) / lb ) ; 
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
% r.volume (m3)
% kinetic energy (J)
% potential energy (J)
% total energy (J)
% vorticity (m2/s)
% enstrophy (J)
% squared normal strain (J)
% squared shear strain (J)
% squared strain (J)
% Okubo-Weiss
% momentum (r.M^3 * r.M / s)a
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

r.u_a;
r.v_a;
r.u;
r.v;
r.u_t;
r.v_t;
r.H;
r.eta;
r.M;
r.N;
r.dx;
r.dy;
r.dt;
rho0;
r.g;
r.K;
r.mask;
r.ke;
r.pe;
r.strechrate_t;
r.shearrate_w;
r.divergence_t;
r.sqstrain_t;
r.sqstrechrate_t;
r.sqshearrate_w;
r.curl_w;
r.enstrophy_t;
r.enstrophy_w;
r.okuboweiss_t;

%Eddy kinetic energy
r.gradux;
r.graduy;
r.gradvx;
r.gradvy;
r.eke;

%r.time-dependent r.properties
r.time;
r.volume;
r.iKe;
r.iPe;
r.vtime;
r.iVort;
r.iEnst;
r.iSqStrech;
r.iSqShear;
r.iSqStrain;
r.iOWeiss;
r.iMomentumU;
r.iMomentumV;
r.iEke;

dA = r.dx * r.dy; %m2

%Computes the average velocity field
r.u_a = ((l - 1) * r.u_a + r.u ) / l;
r.v_a = ((l - 1) * r.v_a + r.v ) / l;

%Eddy kinetic energy
r.gradux = (r.u(2:r.M+1,:) - r.u(1:r.M,:)) / r.dx;
r.gradvy = (r.v(:,2:r.N+1) - r.v(:,1:r.N)) / r.dy;
r.graduy(2:r.M,2:r.N-2) = (r.u(2:r.M,3:r.N-1) - r.u(2:r.M,2:r.N-2)) / r.dy;
r.gradvx(2:r.M-2,2:r.N) = (r.v(3:r.M-1,2:r.N) - r.v(2:r.M-2,2:r.N)) / r.dx;

%Compute diagnostic fields
%Curl
computeCurl_w; % 1 / s

%Strain rate
r.strechrate_t = r.mask .* ((r.u(2:r.M+1,:) - r.u(1:r.M,:)) * r.dy ...
            - (r.v(:,2:r.N+1) - r.v(:,1:r.N)) * r.dx) / dA; % 1 / s
computeShearRate_w % 1 / s

%Growth rate
r.divergence_t = r.mask .* ((r.u(2:r.M+1,1:r.N) - r.u(1:r.M,1:r.N)) * r.dy ...
                + (r.v(1:r.M,2:r.N+1) - r.v(1:r.M,1:r.N)) * r.dx) / dA; % 1/ s

%Potential vorticity
computePotentialVorticity_w

%This is the so-called scalar of Arakawa-Okubo-Weiss.
%Computes the local enstrophy, squared strain rate and Okubo-Weiss
%parameter %Look for Weiss
%La Jolla Institute Report from 1981.
%Physica r.d: Nonlinear 1991.
%Look for Okubo-Weiss parameter/function
%Look for Arakawa 1966 paper.
r.enstrophy_w = .5 * ( r.curl_w.^2); %1 / s2
r.enstrophy_t = interpolFtoT(r.enstrophy_w, r.mask, r.M, r.N);

r.sqshearrate_w = .5 * r.shearrate_w.^2;
r.sqshearrate_t = interpolFtoT(r.sqshearrate_w, r.mask, r.M, r.N);

r.sqstrechrate_t = .5 * r.strechrate_t.^2 .* r.mask;

r.sqstrain_t = r.sqstrechrate_t + r.sqshearrate_t;

r.sqdivergence_t = .5 * r.divergence_t.^2 .* r.mask;

r.okuboweiss_w = (r.sqshearrate_w - r.enstrophy_w);
r.okuboweiss_t = interpolFtoT(r.okuboweiss_w, r.mask, r.M, r.N);
r.okuboweiss_t = r.okuboweiss_t + r.sqstrechrate_t .* r.mask;
%vers�o do Aires do escalar de Okubo-Weiss
r.okuboweiss_t = r.okuboweiss_t  - r.sqdivergence_t;

%%interpolate from r.u and r.v to T cells
r.u_t = .5 * ( r.u(1:r.M,:) + r.u(2:r.M+1,:) );
r.v_t = .5 * ( r.v(:,1:r.N) + r.v(:,2:r.N+1) );

Hnoland = r.H;
Hnoland(find(r.H < -1)) = 0.;

%r.volume
%r.volume(l) = dA * sum(sum(Hnoland)); %m3
%Perturbed r.volume
r.volume(l) = dA * sum(sum(r.eta .* r.mask)); %m3

%Momentum (actually it's only the integrated velocity)
r.iMomentumU(l) = dA * sum(sum(r.u_t .* Hnoland,2)); % kg * r.M / s
r.iMomentumV(l) = dA * sum(sum(r.v_t .* Hnoland,1)); % kg * r.M / s

%Kinetic energy
r.ke = .5 * rho0 * dA * (r.u_t.^2 + r.v_t.^2) .* Hnoland .* r.mask; % kg * r.M * r.M / s2 = J
r.iKe(l) = sum(sum(r.ke)); % kg * r.M * r.M / s2 = J

%Potential energy
r.pe = .5 * rho0 * g * dA .* r.eta.^2 .* r.mask; %kg r.M * r.M / s2 = J
r.iPe(l) = sum(sum(r.pe)); %kg r.M * r.M / s2 = J

%Eddy kinetic energy
r.eke = r.eke + rho0 * dA * r.K * 2 * r.dt * (r.gradux.^2 + r.gradvy.^2 + r.graduy.^2 + r.gradvx.^2) ...
    .* Hnoland .* r.mask;
r.iEke(l) = sum(sum(r.eke));

%r.vorticity
r.iVort(l) = sum(sum(r.curl_w)) * dA; % m2 / s

%r.enstrophy
r.iEnst(l) = sum(sum(r.enstrophy_w)) * dA; % m2 / s2

%r.squared strain
r.iSqShear(l) = sum(sum(r.sqshearrate_w))* dA;
r.iSqStrech(l) = sum(sum(r.sqstrechrate_t)) * dA;
r.iSqStrain(l) = sum(sum(r.sqstrain_t)) * dA;

%r.Okubo-Weiss
r.iOWeiss(l) = sum(sum(r.okuboweiss_t)) * dA; %<--- too large numerical error
%r.iOWeiss(l) = r.iSqStrain(l) - r.iEnst(l);

%r.time
r.vtime(l) = r.time;

diag = [...
        r.volume(l), ...
        r.iKe(l), ...
        r.iPe(l), ...
        r.iVort(l), ...
        r.iEnst(l), ...
        r.iSqStrain(l), ...
        r.iOWeiss(l) ...
        ];

function computeCurl_w
%function computeCurl_w % 1 / s

r.curl_w;
r.curl_t;
r.mask;
r.mask_v;
r.mask_u,
r.u;
r.v;
r.M;
r.N;
r.dy;
r.dx;

dA = r.dx .* r.dy;

%Curl of velocity in the interior...
r.curl_w(2:r.M,2:r.N) = ((r.mask_v(2:r.M,2:r.N) .* r.v(2:r.M,2:r.N) ...
    - r.mask_v(1:r.M-1,2:r.N) .* r.v(1:r.M-1,2:r.N)) .* r.dy ...
    - (r.mask_u(2:r.M,2:r.N) .* r.u(2:r.M,2:r.N) ...
    - r.mask_u(2:r.M,1:r.N-1) .* r.u(2:r.M,1:r.N-1)) .* r.dx) ./ dA; % 1 / s
%... at the corners %Thank you for the conversation with F.L.!
%lower-left SW
r.curl_w(1,1) = ((r.mask_v(1,1) * r.v(1,1)) * r.dy ...
    - (r.mask_u(1,1) * r.u(1,1)) * r.dx)/dA; % 1 / s
%lower-right SE
r.curl_w(1,r.N+1) = ((r.mask_v(1,r.N+1) * r.v(1,r.N+1)) * r.dy ...
    - ( - r.mask_u(1,r.N) * r.u(1,r.N)) * r.dx)/dA; % 1 / s
%top-left NW
r.curl_w(r.M+1,1) = (( - r.mask_v(r.M,1) * r.v(r.M,1)) * r.dy ...
    - (r.mask_u(r.M+1,1) * r.u(r.M+1,1)) * r.dx)/dA; % 1 / s
%top-right NE
r.curl_w(r.M+1,r.N+1) = ((- r.mask_v(r.M,r.N+1) * r.v(r.M,r.N+1)) * r.dy ...
    - (- r.mask_u(r.M+1,r.N) * r.u(r.M+1,r.N)) * r.dx)/dA; % 1 / s
%... at the western boundary
r.curl_w(2:r.M,1) = ((r.mask_v(2:r.M,1) .* r.v(2:r.M,1) - r.mask_v(1:r.M-1,1) .* r.v(1:r.M-1,1)) * r.dy ...
    - (r.mask_u(2:r.M,1) .* r.u(2:r.M,1)) * r.dx)/dA; % 1 / s
%... at the eastern boundary
r.curl_w(2:r.M,r.N+1) = ((r.mask_v(2:r.M,r.N+1) .* r.v(2:r.M,r.N+1) - r.mask_v(1:r.M-1,r.N+1) .* r.v(1:r.M-1,r.N+1)) * r.dy ...
    - (- r.mask_u(2:r.M,r.N) .* r.u(2:r.M,r.N)) * r.dx)/dA; % 1 / s
%... at the southern boundary
r.curl_w(1,2:r.N) = ((r.mask_v(1,2:r.N) .* r.v(1,2:r.N)) * r.dy ...
    - (r.mask_u(1,2:r.N) .* r.u(1,2:r.N) - r.mask_u(1,1:r.N-1) .* r.u(1,1:r.N-1)) * r.dx)/dA; % 1 / s
%... at the northern boundary
r.curl_w(r.M+1,2:r.N) = ((- r.mask_v(r.M,2:r.N) .* r.v(r.M,2:r.N)) * r.dy ...
    - (r.mask_u(r.M+1,2:r.N) .* r.u(r.M+1,2:r.N) - r.mask_u(r.M+1,1:r.N-1) .* r.u(r.M+1,1:r.N-1)) * r.dx)/dA; % 1 / s
% interpolation from the {Z,W,F}-cell to the T-Cell
r.curl_t = interpolFtoT(r.curl_w, r.mask, r.M, r.N);

function computeShearRate_w
%function computeShearRate_w 1 / s
%

r.shearrate_w;
r.shearrate_t;
r.mask;
r.u;
r.v;
r.M;
r.N;
r.dy;
r.dx;

dA = r.dx * r.dy;

%Curl of (-r.u,r.v,0)
%interior of the domain...
r.shearrate_w(2:r.M,2:r.N) = ((r.v(2:r.M,2:r.N) - r.v(1:r.M-1,2:r.N)) * r.dy ...
    + (r.u(2:r.M,2:r.N) - r.u(2:r.M,1:r.N-1)) * r.dx) / dA;
%bottom-left corner (SW)
r.shearrate_w(1,1) = ((r.v(1,1)) * r.dy ...
    + (r.u(1,1)) * r.dx) / dA;
%bottom-right corner (SE)
r.shearrate_w(1,r.N+1) = ((r.v(1,r.N+1)) * r.dy ...
    + (- r.u(1,r.N)) * r.dx) / dA;
%top-left corner (NW)
r.shearrate_w(r.M+1,1) = ((- r.v(r.M,1)) * r.dy ...
    + (r.u(r.M+1,1)) * r.dx) / dA;
%top-right corner (NE)
r.shearrate_w(r.M+1,r.N+1) = ((- r.v(r.M,r.N+1)) * r.dy ...
    + (- r.u(r.M+1,r.N)) * r.dx) / dA;
%western boundary
r.shearrate_w(2:r.M,1) = ((r.v(2:r.M,1) - r.v(1:r.M-1,1)) * r.dy ...
    + (r.u(2:r.M,1)) * r.dx) / dA;
%eastern boundary
r.shearrate_w(2:r.M,r.N+1) = ((r.v(2:r.M,r.N+1) - r.v(1:r.M-1,r.N)) * r.dy ...
    + (- r.u(2:r.M,r.N)) * r.dx) / dA;
%southern boundary
r.shearrate_w(1,2:r.N) = ((r.v(1,2:r.N)) * r.dy ...
    + (r.u(1,2:r.N) - r.u(1,1:r.N-1)) * r.dx) / dA;
%northern boundary
r.shearrate_w(r.M+1,2:r.N) = ((- r.v(r.M,2:r.N)) * r.dy ...
    + (r.u(r.M+1,2:r.N) - r.u(r.M+1,1:r.N-1)) * r.dx) / dA;
% interpolate from F-cell to a T-cell
r.shearrate_t = interpolFtoT(r.shearrate_w, r.mask, r.M, r.N);

function computePotentialVorticity_w
%function computePotentialVorticity_w % 1 / s / r.M
%

r.potvorticity_t;
r.mask;
r.curl_w;
r.potvorticity_w;
r.coriolis;
r.f;
r.H;
r.M;
r.N;

%Adding coriolis
if coriolis
  r.potvorticity_w = r.curl_w + f; % ??? r.curl_w * dA  + f???
end

%potential vorticity
%within the domain...
r.potvorticity_w(2:r.M,2:r.N) = r.potvorticity_w(2:r.M,2:r.N) * 4. ./ ...
    (r.H(2:r.M,2:r.N) + r.H(1:r.M-1,2:r.N) + r.H(1:r.M-1,1:r.N-1) + r.H(2:r.M,1:r.N-1));
%bottom-left corner (SW)
r.potvorticity_w(1,1) = r.potvorticity_w(1,1) / ...
    (r.H(1,1));
%bottom-right corner (SE)
r.potvorticity_w(1,r.N+1) = r.potvorticity_w(1,r.N+1)/ ...
    (r.H(1,r.N));
%top-left corner (NW)
r.potvorticity_w(r.M+1,1) = r.potvorticity_w(r.M+1,1) / ...
    (r.H(r.M,1));
%top-right corner (NE)
r.potvorticity_w(r.M+1,r.N+1) = r.potvorticity_w(r.M+1,r.N+1) / ...
    (r.H(r.M,1));
%western boundary
r.potvorticity_w(2:r.M,1) = r.potvorticity_w(2:r.M,1) * 2. ./ ...
    (r.H(2:r.M,1) + r.H(1:r.M-1,1));
%eastern boundary
r.potvorticity_w(2:r.M,r.N+1) = r.potvorticity_w(2:r.M,r.N+1) * 2. ./ ...
    (r.H(1:r.M-1,r.N) + r.H(2:r.M,r.N));
%southern boundary
r.potvorticity_w(1,2:r.N) = r.potvorticity_w(1,2:r.N) * 2. ./ ...
    (r.H(1,2:r.N) + r.H(1,1:r.N-1));
%northern boundary
r.potvorticity_w(r.M+1,2:r.N) = r.potvorticity_w(r.M+1,2:r.N) * 2. ./ ...
    (r.H(r.M,2:r.N) + r.H(r.M,1:r.N-1));
%interpolate from an F-cell to a T-cell
r.potvorticity_t = interpolFtoT( r.potvorticity_w, r.mask, r.M, r.N);

function tgrid = interpolFtoT(fgrid, r.mask, r.M, r.N)
%% function tgrid = interpolFtoT(fgrid, r.mask, r.M, r.N)
tgrid = .25 * r.mask .* ( ...
              fgrid(1:r.M,1:r.N) ...
            + fgrid(2:r.M+1,1:r.N) ...
            + fgrid(2:r.M+1,2:r.N+1) ...
            + fgrid(1:r.M,2:r.N+1) ...
            );

%--------------------------------------------------------------------------
%SHEL SHallow-water numerical modEL
%Copyright (C) 2006,2009,2010,2011. Guillaume Riflet,
%Instituto Superior T�cnico da Universidade T�cnica de Lisboa.
%--------------------------------------------------------------------------
        