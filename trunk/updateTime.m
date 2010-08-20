function updateTime
%function updateTime

global dt;
global duration;
global outputdt;
global L;
global outputL;

L = getL(duration,dt);
outputL = getL(outputdt, dt);

%% getL (duration_L, dt_L)
function L_L = getL (duration_L, dt_L)
%function L_L = getL (duration_L, dt_L)

    L_L = floor(duration_L / dt_L);