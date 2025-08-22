classdef BathymetryUtilities
    methods (Static)
        function bath = makestep(bath, depth)
            % Creates a stepped bathymetry
            % Example: set half domain to depth
            [M,N] = size(bath);
            bath(:,1:floor(N/2)) = depth;
        end
        function depth = makeisland(depth, x_L, y_L, r)
            % Add island to bathymetry
            [M,N] = size(depth);
            [X,Y] = meshgrid(1:M,1:N);
            mask = ((X-x_L).^2 + (Y-y_L).^2) < r^2;
            depth(mask) = -99; % land value
        end
    end
end
