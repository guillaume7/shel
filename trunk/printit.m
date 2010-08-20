function printit(frame_L, handles, choice)
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
