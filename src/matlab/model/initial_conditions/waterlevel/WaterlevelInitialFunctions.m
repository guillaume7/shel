classdef WaterlevelInitialFunctions
    methods(Static)
        function field = constantField(M, N, value)
            field = value * ones(M, N);
        end
        function field = nullField(M, N)
            field = zeros(M, N);
        end
        function field = gaussianBump(M, N, x0, y0, sx, sy, amplitude, dx, dy)
            [X, Y] = meshgrid((1:N)*dx, (1:M)*dy);
            field = amplitude * exp(-((X-x0).^2/(2*sx^2) + (Y-y0).^2/(2*sy^2)));
        end
        function field = stepField(M, N, value1, value2, stepIndex)
            field = value1 * ones(M, N);
            field(:, stepIndex:end) = value2;
        end
        function field = islandField(M, N, x0, y0, r, landValue, dx, dy)
            [X, Y] = meshgrid((1:N)*dx, (1:M)*dy);
            mask = ((X-x0).^2 + (Y-y0).^2) < r^2;
            field = zeros(M, N);
            field(mask) = landValue;
        end
        function field = pressureDeviation(M, N, value)
            field = value * ones(M, N);
        end
        function field = geostrophicWaterLevel(M, N, eta0, gradient)
            field = eta0 + gradient * repmat((1:N), M, 1);
        end
    end
end
