%--------------------------------------------------------------------------
%     This file is part of SHEL SHallow-water numerical modEL
% 
%     SHEL is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     Foobar is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
%--------------------------------------------------------------------------
function ComputeModel(handles)
%function runmodel_all

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
%
%       T(eta, H, tr, d, f, mask,   U(u, mask_u         V(v, mask_v
%           x, y, wind and bottom      x_u, y_u)           x_v, y_v)
%            stress)
%

%Time parameters and output parameters
global L;
global outputL;
global dt;
global time;
global timetr;
global frame;
global printG;
global film;
global mov;
global myfullfile;
global mydir;
global optprint;

%tracer parameters
global tracer;
global Tr;
global K;

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


%% Write output of initial condition to disk.
if printG    
    if movie
        %Open a new avi file.
        filename = sprintf('%s.avi', [mydir,'/',myfullfile]);
        mov = avifile(filename,'fps', 10,'quality',100,'compression','None');
    end
    printit(frame, handles, optprint);
    frame = frame + 1;
end 

%% Main loop in time
ltr = 0;
for l = 1 : L           
    
    %Update the time
    time = time + dt;
    
    %Compute the time scheme of the momentum+continuity model
    ComputeLeapfrog;
    
    %Numerics to compute the diagnostic variables
    ComputeDiagnostics(l);

    %Compute all the tracers
    if tracer
        if timetr <= time
            [Tr, dttr] = ComputeTracer_FT(Tr, K);
            timetr = time + dttr;
            ltr = ltr + 1;
        end
    end
       
%% Take care of plotting the outputs and write them to disk.
    if ol == outputL
        
        ol=0;
               
        if film
            plotmodel(handles);
            pause(0.1); %s
        end
        if printG
            printit(frame, handles, optprint);
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

%% Take care of writing the final output on disk
if printG
if movie
    mov = close(mov);
end
end

if ~film
    plotmodel(handles);
end

%% Display simulation statistics and performance
cput = toc;
disp(sprintf('Elapsed time: %0.5f s\nTime per iteration: %0.5f s\n', cput, cput/L));
if tracer
    disp(sprintf('Made %0.3d tracer iterations on a total of %0.3d iterations', ltr, l));
end
%--------------------------------------------------------------------------
%SHEL SHallow-water numerical modEL
%Copyright (C) 2006,2009,2010. Guillaume Riflet, Instituto Superior Técnico
%da Universidade Técnica de Lisboa.
%--------------------------------------------------------------------------
    