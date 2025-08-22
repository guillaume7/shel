classdef BottomInitialFunctions
    methods(Static)
        function field = stepField(M, N, value1, value2, stepIndex)
            field = value1 * ones(M, N);
            field(:, stepIndex:end) = value2;
        end
        function field = islandField(M, N, x0, y0, r, landValue, dx, dy)
            [X, Y] = meshgrid((1:N)*dx, (1:M)*dy);
            mask = ((X-x0).^2 + (Y-y0).^2) < r^2;
            field = zeros(M, N); % Default value is 0 (water)
            field(mask) = landValue;
        end
    end
end
