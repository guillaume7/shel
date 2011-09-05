function s = GUI_common_handles
%GUI_COMMON_HANDLES Summary of this function goes here
% Call as:
% s = GUI_common_handles
% y2 = s.GUI_loadsave;
% y3 = s.GUI_initializeglobals;
% ...

s.GUI_loadsave = @loadsave;
s.GUI_initializedglobals = @initializeglobals;
s.GUI_updatestability = @updatestability;
s.GUI_showupvalues = @showupvalues; %must work on it
s.GUI_updatebumpenergy = @updatebumpenergy; %must work on it
s.GUI_UpdateCaracteristicNumbers = @UpdateCaracteristicNumbers; %must work on it
s.GUI_update_options_causal_relationships = @update_options_causal_relationships; %must work on it
s.GUI_plotmodel = @GUI_plotmodel;
s.GUI_printit = @GUI_printit;
s.GUI_setPanels = @setPanels;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pre-callback functions. Local functions that auxiliates the
% menu callbacks.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load and save function
function loadsave(atype, afile)
%function loadsave(type, file)
%
%type - 'load' or 'save'
%Initialize global variables script
%       
%

%Momentum Forces
global coriolis; %use coriolis force? Logical.
global bottom; %use bottom stress? Logical.
global wind; %use wind stress? Logical.
global pressure; %use pressure gradient force? Logical.
global noslip; %use no-slip? Logical.

%Momentum parameters
global rho_air; %air density. Real.
global rho0; %water density. Real.
global Nu; %water turbulent horizontal viscosity. Real.
global g; %gravitational pull. Real.
global f; %Coriolis frequency. Real.
global lb; %Bottom rugosity length. Real.
global uwind; %x-axis wind speed component. Real.
global vwind; %y-axis wind speed component. Real.
global karman; %Von Karman constant. Real.
global gama; %gama coefficient for Asselin-Robert filter. Hidden. Real.

%Time parameters
global duration; %time duration. Real
global dt;      %time step. Real.

%Horizontal Grid parameters
global dx; %x-axis grid step. Real.
global dy; %y-axis grid step. Real.

%View parameters
global outputL;     %number of iterations for output. Hidden. Integer.
global outputdt;    %time interval between graphical outputs. Real.

%Bathymetry
global d0; %Constant depth. Real.
global M; % x-axis grid size. Integer.
global N; % y-axis grid size. Integer
global L; % Characteristic length. Real.
 %Bathymetry file
global loadbathymetry; %bathymetry test. Logical.
global bathymetryfile; %bathymetry file. String.
global d;       %The current bathymetry. Matrix.
global eta0;    % Initial water elevation. Matrix.
 %step bottom
global step;    %Step test. Logical.
global d0_step; %Step depth. Real.
 %island
global tc_isla; %Use island? Logical.
global isla_x0; %island x-position. Real.
global isla_y0; %island y-position. Real.
global isla_R;  %island radius. Real.

%Initial conditions
global tc_taylor; %Use taylor column initialization? Logical.
global tc_geostrophic; %user geostrophic initial fields? Useless. Logical.
 %bump
global tc_bump; %Use gaussian elevation? Logical.
global bump_x0; %gaussian center x-position. Real.
global bump_y0; %gaussian center y-position. Real.
global bump_sx; %gaussian width along x-axis. Real.
global bump_sy; %gaussian width along y-axis. Real.
global bump_d0; %gaussian max height. Real.

%tracer
global tracer; %Use tracer? Logical.
global K; %tracer horizontal turbulent diffusivity.
global TrXo; %tracer center x-position. Real.
global TrYo; %tracer center y-position. Real.
global TrR; %tracer radius. Real.
global Treshold; %tracer min concentration, above which integration is 
                 %made in diagnostics. Real.

%boundary conditions
global radiatelevel; %Use level radiation? Logical. 
global flather; %Use flather radiation? Logical.
global neumann; %Use neumann condition? Logical.
global radiatetan; %Use radiation for tangential velocity? Logical.
global u_closed; %Close the x-axis boundary? Logical.
global v_closed; %Close the y-axis boundary? Logical.

%Panels outputs & properties
global npanels; %Two or four panels? Integer.
global leftprop;  %Property number from property list. For the left panel. Integer.
global rightprop; %Property number from property list. For the right panel. Integer.
 %CLim coef A: [-1eA 1eA]
global CLimEta; % m [-10. 2.] Water level Z-scale limits. Real vector.
global CLimV; % m/s [0. 100.] velocity Z-scale limits. Real vector.
global CLimLeft; % J [-10. 10.] left property z-scale limits. Real vector.
global CLimRight; % J [-10. 10.] right property z-scale limits. Real vector.
 %print to file
global film; %Do we animate the simulation and render it on screen? Logical.
global printG; %Do we export the screen rendering to a file? Logical.
global myfile; %filename to build sub-folder and image filenames. String.
global optprint; %Which export format and which property? Integer.
global movie; %Is it the avi format we want to export? Logical.
global frame; %frame counter of the graphical exports. integer.

disp(nargchk(1, 2, nargin));

switch atype
    case 'load'
        % Kill initial conditions
        clear globals;
        %Eventually change the uiload with
        %the function`*file = uigetfile* then *load(file)*
        uiload;
    case 'save'
        %Eventually change the uisave with
        %the function`*file = uiputfile* then *save(file)*
        switch nargin 
            case 2
                save(afile);
            otherwise
                uisave;
        end
    otherwise
end

%initialize global values
function initializeglobals()
%function initializeglobals()
%
%
global rho0;
global rho_air;
global Nu;
global g;
global f;
global lb;
global karman;
global uwind;
global vwind;
global M;
global N;
global L;
global outputL;     %number of iterations for output
global outputdt;
global dx;
global dy;
global dt;
global duration;
global gama;        %gama coefficient for Asselin-Robert filter
global coriolis;    %use coriolis force?
global bottom;      %use bottom stress?
global wind;        %use wind stress?
global pressure;    %use pressure gradient force?
global noslip;      %use no-slip?

%T-cells
global eta0;

%print to file
global myfile;
global frame;
global movie;
global film;
global printG;
global optprint;

%CLim coef A: 
global CLimEta; % m
global CLimLeft; % J
global CLimRight; % J

%Testcases
global tc_bump;
global tc_geostrophic;
global tc_isla;
global tc_taylor;

%bathymetry
global loadbathymetry;
global bathymetryfile;
global step;
global d0;
global d0_step;

%boundary conditions
global u_closed;
global v_closed;
global u_cyclic;
global v_cyclic;

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
global K;
global TrXo;
global TrYo;
global TrR;
global Treshold;

rho0 = 1000;
rho_air = 1;
Nu = 0;
K = 0;
g = 9.8;
f = 0.;
lb = 0.;
karman = 0.4;
uwind = 0.;
vwind = 0.;
M = 10;
N = 10;
L = 3;
outputL = 1;     %number of iterations for output
outputdt = 2;
dx = 20;
dy = 20;
dt = 10;
duration = 100;
gama = 0.1;        %gama coefficient for Asselin-Robert filter
coriolis = false;    %use coriolis force?
bottom = false;      %use bottom stress?
wind = false;        %use wind stress?
pressure = true;    %use pressure gradient force?
noslip = false;      %use no-slip?

%T-cells
eta0 = 0;

%print to file
myfile = 'null.txt';
frame = true;
movie = false;
film = false;
printG = false;
optprint = false;

%CLim coef A: [-1eA 1eA]
CLimEta = 0.01; % m
CLimLeft = 1.; % J
CLimRight = 1.; % J

%Testcases
tc_bump = false;
tc_geostrophic = false;
tc_isla = false;
tc_taylor = false;

%bathymetry
loadbathymetry = false;
bathymetryfile = '';
step = false;
d0 = 10.;
d0_step = 1.;

%boundary conditions
u_closed = false;
v_closed = false;
u_cyclic = false;
v_cyclic = false;

%bump
bump_x0 = 0;
bump_y0 = 0;
bump_sx = 0;
bump_sy = 0;
bump_d0 = 0.01;

%island
isla_x0 = 0;
isla_y0 = 0;
isla_R = 0;

%tracer
tracer = false;
K = 1.;
TrXo = 0.;
TrYo = 0.;
TrR = 0.;
Treshold = 0.001;

%update stability
function updatestability(handles)
%function updatestability(handles)
%GR : updates the stability status

    global dt;
    global dy;
    global dx;
    global d0;
    global d0_step;
    global g;
    global Nu;
    global tracer;
    global K;

    cfl_U = CFL(dt,dx,dy,max(d0,d0_step),g,Nu);
    
    if tracer
        cfl_T = CFL_T(dt,dx,dy,max(d0,d0_step),g,K);
        set(handles.textTrAC, 'String', num2str( cfl_T(1), '%0.3g'));
        set(handles.textTrB, 'String', num2str( cfl_T(2), '%0.3g'));
        set(handles.textTrUmax, 'String', num2str( cfl_T(3), '%0.3g'));
    end
    
    if cfl_U(1)
        strin = 'Ok';
    else
        strin = 'CFL criterion violation';
    end
     
    %GR : sets the numerical status
    set(handles.statustext, 'String', strin);
        
    %GR : Sets the CFL number
    set(handles.couranttext, 'String', num2str(cfl_U(2),'%0.3g'));
    set(handles.peclettext, 'String', num2str(cfl_U(3),'%0.3g'));

%Shows up the values
function showupvalues(handles)
%function showupvalues(handles)
% --- Is called by resetbutton ----------------

global rho0;
global rho_air;
global Nu;
global g;
global f;
global lb;
global karman;
global uwind;
global vwind;
global M;
global N;
global L;
global outputL;
global dx;
global dy;
global dt;
global coriolis;
global wind;
global bottom;
global pressure;
global noslip;

%Testcases
global tc_bump;
global tc_geostrophic;
global tc_isla;
global tc_taylor;

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

%bath
global loadbathymetry;
global bathymetryfile;
global step;
global d0;
global d0_step;

%boundary conditions
global u_closed;
global v_closed;
global flather;
global radiatelevel;
global radiatetan;
global neumann;

%tracer
global tracer;
global K;
global TrXo;
global TrYo;
global TrR;
global Treshold;

    updatestability(handles);
    updatebumpenergy(handles);
    
    %GR: Sets the checkboxes
    set(handles.bumpbox, 'Value', tc_bump);
    set(handles.geostrophicbox, 'Value', tc_geostrophic);
    set(handles.checkboxIsla, 'Value', tc_isla);
    set(handles.checkboxTaylorCol, 'Value', tc_taylor);
    
    %GR : Sets the the editable parameters
    set(handles.editrho_air, 'String', num2str(rho_air));
    set(handles.editrho0, 'String', num2str(rho0));
    set(handles.editNu, 'String', num2str(Nu));
    set(handles.editg, 'String', num2str(g));
    set(handles.editlatitude, 'String', num2str(getlatitude(f)));
    set(handles.editlb, 'String', num2str(lb));
    set(handles.editkarman, 'String', num2str(karman));
    set(handles.edituwind, 'String', num2str(uwind));
    set(handles.editvwind, 'String', num2str(vwind));
    set(handles.editM, 'String', num2str(M));
    set(handles.editN, 'String', num2str(N));
    set(handles.editduration, 'String', num2str(L*dt));
    set(handles.editoutputdt, 'String', num2str(outputL*dt));
    set(handles.editdx, 'String', num2str(dx));
    set(handles.editdy, 'String', num2str(dy));
    set(handles.editdt, 'String', num2str(dt));
    set(handles.coriolisbox, 'Value', coriolis);
    set(handles.windbox, 'Value', wind);
    set(handles.bottombox, 'Value', bottom);
    set(handles.pressurebox, 'Value', pressure);
    set(handles.wallbox, 'Value', noslip);
    set(handles.editIslaX, 'String', num2str(isla_x0));
    set(handles.editIslaY, 'String', num2str(isla_y0));
    set(handles.editIslaR, 'String', num2str(isla_R));
    set(handles.editBumpX, 'String', num2str(bump_x0));
    set(handles.editBumpY, 'String', num2str(bump_y0));
    set(handles.edit_bump_d0, 'String', num2str(bump_d0));
    set(handles.editBumpSx, 'String', num2str(bump_sx));
    set(handles.editBumpSy, 'String', num2str(bump_sy));
    set(handles.edit_d0, 'String', num2str(d0));
    set(handles.edit_d0step, 'String', num2str(d0_step));
    set(handles.checkbox_step, 'Value', step);
    set(handles.checkboxUClosed, 'Value', u_closed);
    set(handles.checkboxVClosed, 'Value', v_closed);
    set(handles.bathymetrybutton, 'Value', loadbathymetry);
    set(handles.editbathymetry, 'String', bathymetryfile);
    set(handles.checkboxTracer, 'Value', tracer);
    set(handles.editKtr, 'String', num2str(K));
    set(handles.editTrXo, 'String', num2str(TrXo));
    set(handles.editTrYo, 'String', num2str(TrYo));
    set(handles.editTrR, 'String', num2str(TrR));
    set(handles.editTrThreshold, 'String', num2str(Treshold));    
    set(handles.popupmenuLevelOB,'Value', ...
        radiatelevel * 2 + ~radiatelevel * 1);
    set(handles.popupmenuUOB,'Value', ...
        flather * 3 + neumann * 2 + ~(flather || neumann) * 1);
    set(handles.popupmenuVOB,'Value', ...
        radiatetan * 2 + ~radiatetan * 1)
   
    UpdateCaracteristicNumbers(handles);
    
    update_options_causal_relationships(handles);    

%function lat = getlatitude(f) in degrees
function lat_L = getlatitude(f_L)
%function lat = getlatitude(f) in degrees

    omega_L = .00007272; %rad/s Mellor Intro to Phys Oceanography
    lat_L = 180. / pi * asin( .5 * f_L / omega_L  );
    
%cfl criterion
function cfl_T = CFL_T(dt, dx, dy, H, g, K)
    %Pt+dt = A Pi+1 + B Pi + C Pi-1
    % 0 <= A,B,C <= 1
    %CFL computed for an upwind tracer

    %Estimated max flow velocity, based on the wave
    %wake velocity of the flow for ripples at the limit
    %of the hydrostatic approximation.
    h = H / 100; %Wave height at the limit of the hydrostatic appr.
    v_est =  h / 2 * sqrt(g / H);

    Cr = v_est * dt ./ min([dx dy]);
    Dif = K * dt ./ min([dx dy]).^2;
    AC = Cr + Dif;
    B = 1 - 4 * AC; %Instabilities (UP or CD)
    peclet = Cr / 2 / Dif; %Negative oscillations (CD)
    Umax = .25 * dx / dt - K / dx; %Max vel (UP)
    
    cfl_T = [B peclet Umax];

function cfl_L = CFL(dt_L, dx_L, dy_L, H_L, g_L, nu)
%function bool_L = CFL(dt_L, dx_L, dy_L, H_L, g_L, nu)
%
%CFL criterion enumerated in Kantha-Clayson p.282 ~wrong
%for the shallow water equations solver.

 c = sqrt( g_L * max(max(H_L)) );

 transportivity = 2 * c * dt_L * ( dx_L^(-1) + dy_L^(-1) );
 
 %Let us estimate v_max based on the flow speed in the wake of ripples
 %travelling across the domain instead
 v_max =  H_L / 100 * sqrt(g_L / H_L);
 %v_max = max([dx_L dy_L]) / dt_L;
 
 positivity = 0.5 * v_max * max([dx_L dy_L]) ./ nu;
 
 if transportivity < 1 && positivity < 1
     cfl_L = [true, transportivity, positivity];
 else
     cfl_L = [false, transportivity, positivity];
 end

 function updatebumpenergy(handles)
    global bump_d0;
    global bump_sx;
    global bump_sy;
    global rho0;
    global g;
    global d0;
    set(handles.bumpenergytext, 'String', ...
        num2str(.25 * rho0 * g * bump_d0^2 * bump_sx * bump_sy * pi,'%0.3g'));
    set(handles.bumpUotext, 'String', ...
        num2str(.5 * bump_d0 * sqrt(g/d0),'%0.1e'));

function UpdateCaracteristicNumbers(handles)
%This function updates the Reynolds and the Rossby numbers

    %read these variables
    global f;
    global g; %does not has an update function
    global Nu;
    global d0;
    global d0_bump;
    global dx;
    global N;
    global M;
    global coriolis;

    %update these variables
    global Ro;
    global Re;
    global Celerity;
    global Uo;
    global RoRadius;
    global CaracTime;
    global HLength;
    global VLength;
    
    %Compute the values
    VLength = d0;
    HLength = max(N * dx, M * dx);
    Celerity = sqrt(g * VLength);
   % Uo = .5 * d0_bump * sqrt( g / (d0_bump + d0));
    Uo = 0.1; %m/s
    %Ro = sqrt(g/VLength)/f;
    Ro = Uo / (f * HLength);
    Re = Celerity * VLength / Nu;
    RoRadius = Celerity / f;
    if coriolis
        CaracTime = max( max( sqrt(VLength/g), 1/f ), HLength/Celerity );
    else
        CaracTime = max( sqrt(VLength/g), HLength/Celerity);
    end
    
    %Write the new values
    set(handles.ReynoldsText, 'String', num2str(Re,'%0.3g'));
    set(handles.CelerityText, 'String', num2str(Celerity,'%0.3g'));
    set(handles.CaracTimeText, 'String', num2str(CaracTime,'%0.3g'));
    set(handles.HSizeText, 'String', num2str(HLength,'%0.3g'));
    set(handles.VSizeText, 'String', num2str(VLength,'%0.3g'));
    if coriolis
        set(handles.RossbyText, 'String', num2str(Ro,'%0.3g'));
        set(handles.RossbyRadiusText, 'String', num2str(RoRadius,'%0.3g'));
    else
        set(handles.RossbyText, 'String', '-');
        set(handles.RossbyRadiusText, 'String', '-');
    end

%-This function will update all the causal relationships!        
function update_options_causal_relationships(handles)
    
    global u_closed;
    global v_closed;
    global loadbathymetry;
    global coriolis;
    global pressure;
    global wind;
    global bottom;
    global tc_bump;
    global tc_isla;
    global tc_taylor;
    global step;
    global tracer;
    
    %u-cyclic or not? UNCOMMENT THIS IF YOU WANT TO KEEP DEVELOPPING
    %CYCLIC BOUNDARIES
%    if (v_closed)
%        de_activate_v_cyclic(handles, 'off')
%    else
%        de_activate_v_cyclic(handles, 'on')
%    end
    
    %v-cyclic or not? UNCOMMENT THIS IF YOU WANT TO KEEP DEVELOPPING
    %CYCLIC BOUNDARIES
%    if (u_closed)
%        de_activate_u_cyclic(handles, 'off')
%    else
%        de_activate_u_cyclic(handles, 'on')
%    end

    %Do we need OBC?
    if (u_closed && v_closed)
        de_activate_obc(handles, 'off')
    else
        de_activate_obc(handles, 'on')
    end

    %manual bathymetry?
    if (loadbathymetry)
        de_activate_bathymetry(handles, 1);
    else
        de_activate_bathymetry(handles, 0);
    end
    
    %Coriolis?
    if coriolis
        de_activate_coriolis(handles, 'on')
    else
        de_activate_coriolis(handles, 'off')
    end

    %Pressure gradient force?
    if pressure
        de_activate_pressure(handles,'on')
    else
        de_activate_pressure(handles,'off')
    end
    
    %wind?
    if wind
        de_activate_wind(handles,'on');
    else
        de_activate_wind(handles,'off');
    end
    
    %bottom?
    if bottom
        de_activate_bottom(handles,'on');
    else
        de_activate_bottom(handles,'off');
    end
    
    %bump or taylor?
    if (tc_bump || tc_taylor)
        de_activate_bump(handles, 'on')
    else
        de_activate_bump(handles, 'off')
    end
    
    %island?
    if (tc_isla)
        de_activate_isla(handles, 'on');
    else
        de_activate_isla(handles, 'off');
    end
    
    %step?
    if step
        de_activate_step(handles, 'on');
    else
        de_activate_step(handles, 'off');
    end
    
    %tracer?
    if tracer
        de_activate_tracer(handles, 'on');
    else
        de_activate_tracer(handles, 'off');
    end
    
function setPanels(npanels, handles)
%function setPanels(npanels, handles)
%2 --> 2 panels
%4 --> 4 panels
    if npanels == 2
        pan2 = true;
    else
        pan2 = false;
    end
%    set(handles.radiobutton2pan, 'Value', pan2);
%    set(handles.radiobutton4pan, 'Value', ~pan2);

%% model output plotting
function GUI_plotmodel(handles)
%function plotmodel(handles)

    global CLimEta;
    global CLimLeft;
    global CLimRight;
    global leftprop;
    global rightprop;
    global npanels;
    global time;
    
    %Plot the time
    strin = sprintf('time = %0.3f s', time);
    set(handles.timetext,'String',strin);

    %Plot the waterlevel
    plotproperty(2,CLimEta,handles.levelaxes);
    
    %Plot the flow vector field AND superimpose the bathymetry contours
    plotproperty(4,0,handles.velocityaxes);
    
    if npanels == 4
      
        plotproperty(leftprop,CLimLeft,handles.Leftaxes);
        
        plotproperty(rightprop,CLimRight,handles.Rightaxes);
        
    end

function plotproperty(prop,CLim,Hn)
%function plotproperty(prop,CLim,Hn)
%
%'prop' available list of values:
%  1 - bathymetry
%  2 - eta
%  3 - uv
%  4 - Hu, Hv
%  5 - T mask
%  6 - U mask
%  7 - V mask
%  8 - W mask
%  9 - u
% 10 - v
% 11 - velocity modulus
% 12 - ke
% 13 - pe
% 14 - e
% 15 - eke
% 16 - curl (w)
% 17 - enstrophy (w^2)
% 18 - squared shear (sigma_s^2)
% 19 - squared strech (sigma_n^2)
% 20 - squared strain (sigma^2)
% 21 - Okubo-Weiss (W)
% 22 - global volume
% 23 - global UV-momentum
% 24 - global e-ke-pe-eke
% 25 - global vorticity (w)
% 26 - global enstrophy (w^2)
% 27 - global squared shear
% 28 - global squared strech
% 29 - global squared strain
% 30 - global OW (W, w^2, sigma^2)
% 31 - YZ Section level
% 32 - Okubo-Weiss (contours)
% 33 - Tracer

    global M;
    global N;
    global bump_d0;
    global x;
    global y;
    global x_u;
    global y_u;
    global x_v;
    global y_v;
    global x_w;
    global y_w;
    global eta;
    global H;
    global Tr;
    global u_t;
    global v_t;
    global curl_t;
%    global potvorticity_t;
    global okuboweiss_t;
    global sqshearrate_t;
    global sqstrechrate_t;
    global sqstrain_t;
    global enstrophy_t;
    global masknan;
    global mask;
    global mask_u;
    global mask_v;
    global mask_z;
    global d;
    global ke;
    global pe;
    global eke;
    global vtime;
    global iVort;
    global iEnst;
    global iSqStrech;
    global iSqShear;
    global iSqStrain;
    global iOWeiss;
    global iMomentumU;
    global iMomentumV;
    global iKe;
    global iPe;
    global iEke;
    global volume;
%    global lb; %bottom roughness
%    global d0; %mean depth
    
    global fontsize;
    global titfontsize;
    global legfontsize;
    
    fontsize = 14;
    titfontsize = fontsize + 4;
    legfontsize = fontsize - 2;
    
    %Store the view (useful for the surf fields)
    [az,el] = view(Hn);
    
        switch prop

            case 1 %Bathymetry (d)
                pcolor(Hn,x,y,d);
                shading(Hn,'interp');
                title(Hn,'Bathymetry','FontSize',titfontsize);
                SetXY(Hn);            
                colorbar('peer',Hn,'FontSize',fontsize);
                            
            case 2 %Eta
                surf(Hn,x(2:M-1,2:N-1),y(2:M-1,2:N-1), masknan(2:M-1,2:N-1) .* eta(2:M-1,2:N-1));
                %title(Hn,'Water Level');
                SetXY(Hn);
                SetZ(CLim,'Water Level, \it \eta \rm (m)',Hn, az, el);                
                %SetLighting(S);

            case 3 %uv
                nn = 30; %Maximum of 30 arrows per axis
                s = ceil(M / nn);
                t = ceil(N / nn);
                quiver(Hn,x(2:s:M-1,2:t:N-1),y(2:s:M-1,2:t:N-1), ...
                    masknan(2:s:M-1,2:t:N-1) .* u_t(2:s:M-1,2:t:N-1), ...
                    masknan(2:s:M-1,2:t:N-1) .* v_t(2:s:M-1,2:t:N-1));
                if min(min(d(2:M-1,2:N-1))) ~= max(max(d(2:M-1,2:N-1)))
                    hold(Hn,'on');
                    contour(Hn,x(2:M-1,2:t:N-1),y(2:M-1,2:t:N-1),d(2:M-1,2:t:N-1));
                    hold(Hn,'off');
                end
                title(Hn,'Velocity Field, (\itu\rm, \itv\rm)','FontSize',titfontsize);
                %legend(Hn,'velocity vector (m/s)');
                SetXY(Hn);

            case 4 %Hu, Hv
                nn = 30; %Maximum of 30 arrows per axis
                s = ceil(M / nn);
                t = ceil(N / nn);
                quiver(Hn,x(2:s:M-1,2:t:N-1),y(2:s:M-1,2:t:N-1), ...
                    masknan(2:s:M-1,2:t:N-1) .* H(2:s:M-1,2:t:N-1) .* u_t(2:s:M-1,2:t:N-1), ...
                    masknan(2:s:M-1,2:t:N-1) .* H(2:s:M-1,2:t:N-1) .* v_t(2:s:M-1,2:t:N-1));
                if min(min(d(2:M-1,2:N-1))) ~= max(max(d(2:M-1,2:N-1)))
                    hold(Hn,'on');
                    contour(Hn,x(2:M-1,2:t:N-1),y(2:M-1,2:t:N-1),d(2:M-1,2:t:N-1));
                    hold(Hn,'off');
                end
                title(Hn,'Flow Field, (\itHu\rm, \itHv\rm)','FontSize',titfontsize);
                %legend(Hn,'velocity vector (m/s)');
                SetXY(Hn);
            
            case 5 %Land mask
                pcolor(Hn,x,y,-(mask-1));
                shading(Hn,'interp');
                title(Hn,'Land Mask','FontSize',titfontsize);
                SetXY(Hn);            
            
            case 6 %U mask
                pcolor(Hn,x_u,y_u,-(mask_u-1));
                shading(Hn,'interp');
                title(Hn,'U Mask','FontSize',titfontsize);
                SetXY(Hn);            
            
            case 7 %V mask
                pcolor(Hn,x_v,y_v,-(mask_v-1));
                shading(Hn,'interp');
                title(Hn,'V Mask','FontSize',titfontsize);
                SetXY(Hn);            
            
            case 8 %W mask
                pcolor(Hn,x_w,y_w,-(mask_z-1));
                shading(Hn,'interp');
                title(Hn,'W Mask','FontSize',titfontsize);
                SetXY(Hn);            
          
            case 9 %u
                surf(Hn,x(2:M-1,2:N-1),y(2:M-1,2:N-1), masknan(2:M-1,2:N-1) .* u_t(2:M-1,2:N-1));
                SetXY(Hn);
                SetZ(CLim,'u Velocity, \itu\rm (m s^{-1})',Hn, az, el);

            case 10 %v
                surf(Hn,x(2:M-1,2:N-1),y(2:M-1,2:N-1), masknan(2:M-1,2:N-1) .* v_t(2:M-1,2:N-1));
                SetXY(Hn);
                SetZ(CLim,'v Velocity, \itv\rm (m s^{-1})',Hn, az, el);

            case 11 %velocity modulus
                surf(Hn,x(2:M-1,2:N-1),y(2:M-1,2:N-1), masknan(2:M-1,2:N-1) .* sqrt( u_t(2:M-1,2:N-1).^2 + v_t(2:M-1,2:N-1).^2 ));
                SetXY(Hn);              
                SetZ(CLim,'Velocity Modulus, (m s^{-1})',Hn, az, el);
                
            case 12 %ke
                surf(Hn,x(2:M-1,2:N-1),y(2:M-1,2:N-1), masknan(2:M-1,2:N-1) .* ke(2:M-1,2:N-1));
                %title(Hn,'Kinetic energy');
                SetXY(Hn);
                SetZ(CLim,'Kinetic Energy Per Grid-Cell, \itK\rm (J)',Hn,az,el);

            case 13 %pe
                surf(Hn,x(2:M-1,2:N-1),y(2:M-1,2:N-1), masknan(2:M-1,2:N-1) .* pe(2:M-1,2:N-1));
                %title(Hn,'Available potential energy');
                SetXY(Hn);
                SetZ(CLim,'Potential Energy Per Grid-Cell,  \itP\rm (J)',Hn,az,el);

            case 14 %e
                surf(Hn,x(2:M-1,2:N-1),y(2:M-1,2:N-1), masknan(2:M-1,2:N-1) .* ( pe(2:M-1,2:N-1) + ke(2:M-1,2:N-1) ));
                %title(Hn,'Energy');
                SetXY(Hn);              
                SetZ(CLim,'Total Energy Per Grid-Cell,  \itTE\rm (J)',Hn,az,el);
            
            case 15 %eke
                surf(Hn,x(2:M-1,2:N-1),y(2:M-1,2:N-1), masknan(2:M-1,2:N-1) .* eke(2:M-1,2:N-1) );
                %title(Hn,'Energy');
                SetXY(Hn);              
                SetZ(CLim,'Total Eddy Kinetic Energy Per Grid-Cell,  \itEKE\rm (J)',Hn,az,el);
            
            case 16 %S - Curl
                surf(Hn,x(2:M-1,2:N-1),y(2:M-1,2:N-1), masknan(2:M-1,2:N-1) .* curl_t(2:M-1,2:N-1));
                %title(Hn,'Curl');
                SetXY(Hn);
                SetZ(CLim,'Velocity Curl,  \omega (s^{-1})',Hn,az,el);
            
            case 17 % Enstrophy
                surf(Hn,x(2:M-1,2:N-1),y(2:M-1,2:N-1), masknan(2:M-1,2:N-1) .* enstrophy_t(2:M-1,2:N-1));
                %title(Hn,'Enstrophy');
                SetXY(Hn);
                SetZ(CLim,'Enstrophy, \frac{1}{2}\omega^2 (s^{-2})',Hn,az,el);
                
            case 18 % Squared shear stress
                surf(Hn,x(2:M-1,2:N-1),y(2:M-1,2:N-1), masknan(2:M-1,2:N-1) .* sqshearrate_t(2:M-1,2:N-1));
                %title(Hn,'Squared shear stress');
                SetXY(Hn);
                SetZ(CLim,'Squared Shear Strain Stress,  \frac{1}{2}\sigma_s^2 (s^{-2})',Hn,az,el);
                
            case 19 % Squared strech stress
                surf(Hn,x(2:M-1,2:N-1),y(2:M-1,2:N-1), masknan(2:M-1,2:N-1) .* sqstrechrate_t(2:M-1,2:N-1));
                %title(Hn,'Squared strech stress');
                SetXY(Hn);
                SetZ(CLim,'Squared Normal Strain Stress,  \frac{1}{2}\sigma_n^2 (s^{-2})',Hn,az,el);
                
            case 20 % Squared strain stress
                surf(Hn,x(2:M-1,2:N-1),y(2:M-1,2:N-1), masknan(2:M-1,2:N-1) .* sqstrain_t(2:M-1,2:N-1));
                title(Hn,'Squared Strain Stress');
                SetXY(Hn);
                SetZ(CLim,'Squared Strain Stress,  \frac{1}{2}\sigma^2 (s^{-2})',Hn,az,el);
                
            case 21 % Scalar of Okubo-Weiss (1/s2)
                surf(Hn,x(2:M-1,2:N-1),y(2:M-1,2:N-1), masknan(2:M-1,2:N-1) .* okuboweiss_t(2:M-1,2:N-1));
                %title(Hn,'Okubo-Weiss');
                SetXY(Hn);
                SetZ(CLim,'Okubo-Weiss,  \itW\rm (s^{-2})',Hn,az,el);

            case 22 %global volume variation(m3)
                plot(Hn,vtime, volume);
 %               axis(Hn,[min(vtime) max(vtime) min(volume)*.999999 max(volume)*1.000001]);
                xlabel(Hn,'Time,  \itt\rm (s)', 'FontSize', fontsize);
                ylabel(Hn,'Volume,  \itV\rm (m^3)', 'FontSize', fontsize);
                title(Hn,'Volume deviation', 'FontSize', titfontsize);
                set(Hn, 'FontSize', fontsize);
                
            case 23 %U and V momentum (actually it's velocity only)
                plot(Hn,vtime,iMomentumU, '-r', 'linewidth',1.5);
                hold(Hn, 'on');
                plot(Hn,vtime,iMomentumV, '-.k', 'linewidth',1.3);
                hold(Hn, 'off');
                xlabel(Hn, 'Time, \itt\rm (s)', 'FontSize', fontsize);
                ylabel(Hn,'Integrated Velocity, (m^4 s^{-1})', 'FontSize', fontsize);
                title(Hn,'U and V Integrated Velocity Evolution', 'FontSize', fontsize);
                set(Hn, 'FontSize', fontsize);
                l=legend(Hn,'\it U \rm','\it V \rm');
                set(l,'FontSize',legfontsize);

            case 24 %e-ke-pe-eke
                plot(Hn,vtime, iKe + iPe, '-k', 'linewidth',1.5); hold(Hn, 'on');
                plot(Hn,vtime, iKe, '-r', 'linewidth',1.3);
                plot(Hn,vtime, iPe, '--b', 'linewidth',1.3);
                plot(Hn,vtime, iEke + iKe + iPe, '-.k', 'linewidth',1.1); hold(Hn, 'off');
                xlabel(Hn,'Time,  \itt\rm (s)', 'FontSize', fontsize);
                ylabel(Hn,'Energy,  \itTE, K, P\rm (J)', 'FontSize', fontsize);
                title(Hn,'Eddy, Kinetic, Potential and Total Energy Evolution', 'FontSize', fontsize);
                set(Hn, 'FontSize', fontsize);
                l=legend(Hn,'\itTE\rm','\itK\rm','\itP\rm','\itEKE\rm');
                set(l,'FontSize',legfontsize);
            
            case 25 %integrated vorticity
                if( max(iVort) < 1e-25 )
                    a = 1e25;
                else
                    a = 1;
                end
                plot(Hn, vtime, iVort);
                %axis(Hn, [min(vtime) max(vtime) min(iVort) * a * 0.9999999 max(iVort) * a * 1.0000001]);
                xlabel(Hn, 'Time,  \itt\rm (s)', 'FontSize', fontsize);
                ylabel(Hn, 'Integrated Vorticity (m^2 s^{-1})', 'FontSize', fontsize);
                title(Hn, 'Integrated Vorticity Evolution', 'FontSize', fontsize);
                set(Hn, 'FontSize', fontsize);
            
            case 26 %iEnstrophy (J)
                plot(Hn,vtime, iEnst);
                xlabel(Hn,'Time,  \itt\rm (s)', 'FontSize', fontsize);
                ylabel(Hn,'Integrated Enstrophy (m^2 s^{-2})', 'FontSize', fontsize);
                title(Hn,'Enstrophy Evolution', 'FontSize', fontsize);
                set(Hn, 'FontSize', fontsize);
            
            case 27 %iSqShear
                plot(Hn,vtime, iSqShear);
                xlabel(Hn,'Time,  \itt\rm (s)', 'FontSize', fontsize);
                ylabel(Hn,'Integrated Squared Shear Strain Stress (m^2 s^{-2})', 'FontSize', fontsize);
                title(Hn,'Integrated Squared Shear Strain Stress', 'FontSize', fontsize);
                set(Hn, 'FontSize', fontsize);
            
            case 28 %iSqStrech
                plot(Hn,vtime, iSqStrech);
                xlabel(Hn,'Time,  \itt\rm (s)', 'FontSize', fontsize);
                ylabel(Hn,'Integrated Squared Normal Strain Stress (m^2 s^{-2})', 'FontSize', fontsize);
                title(Hn,'Integrated Squared Normal Strain Stress', 'FontSize', fontsize);
                set(Hn, 'FontSize', fontsize);
            
            case 29 %iSqStrain (J)
                plot(Hn,vtime, iSqStrain);
                xlabel(Hn,'Time,  \itt\rm (s)', 'FontSize', fontsize);
                ylabel(Hn,'Integrated Squared Strain Stress (m^2 s^{-2})', 'FontSize', fontsize);
                title(Hn,'Integrated Squared Strain Stress Evolution', 'FontSize', fontsize);
                set(Hn, 'FontSize', fontsize);
            
            case 30 % global OW
                plot(Hn,vtime, -iEnst, '--b', 'LineWidth',1.3); hold(Hn,'on'); 
                plot(Hn,vtime, iSqStrain, '-r', 'LineWidth',1.3); hold(Hn,'on'); 
                plot(Hn,vtime, iOWeiss, '-.k', 'LineWidth',1.5); hold(Hn,'off');
                xlabel(Hn,'Time,  \itt\rm (s)', 'FontSize', fontsize);
                ylabel(Hn,'-Enstrophy, Squared Strain, Okubo-Weiss (m^2 s^{-2})', 'FontSize', fontsize);
                set(Hn, 'FontSize', fontsize);
                l=legend(Hn,'-Enstrophy','Squared strain','Okubo-Weiss');
                set(l,'FontSize',legfontsize);
                title(Hn,'-Enstrophy, Squared Strain and Okubo-Weiss', 'FontSize', fontsize);

            case 31 %Waterlevel section
                %surf(Hn,x(2:M-1,2:N-1),y(2:M-1,2:N-1), masknan(2:M-1,2:N-1) .* eta(2:M-1,2:N-1));
                plot(Hn,x(2:M-1,round(N/2)), masknan(2:M-1,round(N/2)) .* eta(2:M-1,round(N/2)) );
                axis(Hn, [min(x(2:M-1,round(N/2))) max(x(2:M-1,round(N/2))) -bump_d0*1.05 bump_d0*1.05]);
                xlabel(Hn,'Position Along X-Axis,  \itx\rm (m)', 'FontSize', fontsize);
                ylabel(Hn,'Water Level,  \it\eta\rm (m)', 'FontSize', fontsize);
                title(Hn,'Water Level Cross Section Evolution', 'FontSize', titfontsize);
                set(Hn, 'FontSize', fontsize);

            case 32 %Okubo-Weiss (contours)
                zlines=[.5 4 8]*1e-21;
                contour(Hn, x,y, masknan .* okuboweiss_t,zlines);
                hc = get(Hn,'Children');
                set(hc,'color','k');
                for i=1:length(hc)
                    hcd = get(hc(i),'UserData');
                    if (hcd>0)
                        set(hc(i),'LineWidth',1.5);
                        set(hc(i),'LineStyle',':');
                    elseif (hcd==0)
                        set(hc(i),'LineWidth',3);
                    else %if (hcd<0)
                        set(hc(i),'LineWidth',2);
                    end
                end
                SetXY(Hn);
                set(Hn,'YLim', [min(min(y)) max(max(y))] );
                set(Hn,'XLim', [min(min(x)) max(max(x))] );
                title(Hn,'Okubo-Weiss parameter');

            case 33 % Tracer plot
                pcolor(Hn,x,y,Tr);
                shading(Hn,'interp');
                title(Hn,'Tracer concentration','FontSize',titfontsize);
                SetXY(Hn);
                colorbar('peer',Hn,'FontSize',fontsize);
                
            otherwise

        end
                
function SetXY(Hn)
%function SetXY - sets the X and Y scales

    global x;
    global y;
    global fontsize;

    %Fonsize of tickLabels
    set(Hn, 'FontSize', fontsize);
    
    %Y
    set(Hn,'YLimMode','manual');
    set(Hn,'YLim', [min(min(y)) max(max(y))] );
    ylabel(Hn,'Position,\it y \rm (m)', 'FontSize', fontsize);

    %X
    set(Hn,'XLimMode','manual');
    set(Hn,'XLim', [min(min(x)) max(max(x))] );
    xlabel(Hn,'Position, \it x \rm (m)', 'FontSize', fontsize);

function SetZ(CLim, label, Hn, az, el)
%function SetZ - sets the Z scale

    global fontsize;
    
    %Color
    set(Hn,'CLimMode','manual');
    set(Hn,'CLim', [-10^(CLim) 10^(CLim)] );
    colormap(Hn,winter);
    colorbar('peer',Hn,'FontSize',fontsize);

    %View
    %set(Hn,'View',[az_level el_level],'Color','k','proj','p');
    view(Hn,[az,el]);
    set(Hn,'Color','k');    
    
    %Light
%GREAT LIGHTING PARAMETERS!    
%lightangle(L, 60, 90);   
%set(gcf,'Renderer','zbuffer')    
    L=light('Parent',Hn);
    lightangle(L, 60, 90);   
    set(gcf,'Renderer','zbuffer')    
    
    %material('shiny');
    set(Hn ...
        ,'DefaultSurfaceAmbientStrength', 0.3 ...
        ,'DefaultSurfaceDiffuseStrength', 0.6 ...
        ,'DefaultSurfaceSpecularStrength', 0.9 ...
        ,'DefaultSurfaceSpecularExponent', 20 ...
        ,'DefaultSurfaceSpecularColorReflectance', 1.0 ...
    );
    shading(Hn, 'interp');
    lighting(Hn, 'phong');

    %Z
    set(Hn,'ZLimMode','manual');
    set(Hn,'ZLim', [-10^(CLim) 10^(CLim)] );
    zlabel(Hn,label, 'FontSize', fontsize);   

function SetLighting(Hn)
%function SetLighting(Hn)

   %alpha(Hn, 0.9);
   
%% model output printing   
function GUI_printit(frame_L, handles, choice)
%function printit
%
%  1 - png_all
%  2 - png_level
%  3 - png_vel
%  4 - png_left
%  5 - png_right
%  6 - eps_all
%  7 - eps_level
%  8 - eps_vel
%  9 - eps_left
% 10 - eps_right
% 11 - avi_all
% 12 - avi_level
% 13 - avi_vel
% 14 - avi_left
% 15 - avi_right
%    
    global mydir;
    global myfullfile;
    global mov;

%    filename = cat(2,cat(2,myfile,'_'),num2str(frame_L));
    filename = sprintf('%s-%0.3d', myfullfile, frame_L);
    move = false;

    %Define the type of output
    switch choice

        case {1, 2, 3, 4, 5}
            typeo = '-dpng';
            dpires = '-r600'; %600 dpi!!
        case {6, 7, 8, 9, 10}
            typeo = '-depsc';
            dpires = '-r300';
        case {11, 12, 13, 14, 15}
            move = true;

    end

    %Define what figure to print
    switch choice

        case {1, 6, 11} % all
            currfig = handles.figure1;

        case {2, 7, 12} % level
            curraxes = handles.levelaxes; 

        case {3, 8, 13} % velocity
            curraxes = handles.velocityaxes; 

        case {4, 9, 14} % left
            curraxes = handles.Leftaxes;

        case {5, 10, 15} % right
            curraxes = handles.Rightaxes; 

    end

    %If a new figure is required.
    switch choice

        case {1, 6, 11} % all

        otherwise % level, velocity, left and right
            currfig = figure;
            currlegend = legend(curraxes);            
            curraxes = copyobj(curraxes, currfig);
%            set(curraxes,'ActivePositionProperty', 'OuterPosition');
            set(curraxes,'outerposition', [3.5 2.5 150 40]);
            colormap(curraxes, winter);
            if length(currlegend)==1
                legend(curraxes, get(currlegend,'String'));
            end
    end

    %In some cases a new figure was required.
    if move
       mov = addframe(mov, getframe(currfig));
    else
       print(currfig, typeo, dpires, [mydir,'/',filename]);
    end
    
    %If a new figure was created, then close it.
    switch choice

        case {1, 6, 11} % all

        otherwise % level, velocity, left and right
            close(currfig);

    end
    
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END OF Pre-callback functions.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

