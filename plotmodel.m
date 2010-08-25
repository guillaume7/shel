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
function plotmodel(handles)
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


%% function plotproperty(prop,CLim,Handle)
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
                
%% function SetXY - sets the X and Y scales
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

%% function SetZ - sets the Z scale
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

%% function SetLighting(Hn)
function SetLighting(Hn)
%function SetLighting(Hn)

   %alpha(Hn, 0.9);
            

%--------------------------------------------------------------------------
%SHEL SHallow-water numerical modEL
%Copyright (C) 2006,2009,2010. Guillaume Riflet, Instituto Superior Técnico
%da Universidade Técnica de Lisboa.
%--------------------------------------------------------------------------
