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
%--------------------------------------------------------------------------%After a breakpoint in the plotmodel.m
figure;
gcf;gca;
ff=20;
ll=54;
zlines=[.5 1 2]*1e-21;
[c,h]=contour(x(ff:ll,ff:ll),y(ff:ll,ff:ll), masknan(ff:ll,ff:ll) .* okuboweiss_t(ff:ll,ff:ll),zlines);
hc = get(h,'Children');
set(h,'color','k');
for i=1:length(hc)
    hcd = get(hc(i),'UserData');
    if (hcd>0);
        set(hc(i),'LineWidth',1.5);
        set(hc(i),'LineStyle',':');
    end
    if (hcd==0);
        set(hc(i),'LineWidth',3);
    end
    if (hcd<0);
        set(hc(i),'LineWidth',2);
    end
end
hold on; gcf;gca;
%pch=pcolor(x(ff:ll,ff:ll),y(ff:ll,ff:ll), masknan(ff:ll,ff:ll) .* okuboweiss_t(ff:ll,ff:ll));
gcf;gca;
shading interp;
%colorbar; colormap summer;
SetXY(gca);
set(gca,'YLim', [min(min(y(ff:ll,ff:ll))) max(max(y(ff:ll,ff:ll)))] );
set(gca,'XLim', [min(min(x(ff:ll,ff:ll))) max(max(x(ff:ll,ff:ll)))] );
title('Okubo-Weiss parameter');
hold off;

%--------------------------------------------------------------------------
%SHEL SHallow-water numerical modEL
%Copyright (C) 2006,2009,2010. Guillaume Riflet, Instituto Superior Técnico
%da Universidade Técnica de Lisboa.
%--------------------------------------------------------------------------
