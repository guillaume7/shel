%After a breakpoint in the plotmodel.m
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
