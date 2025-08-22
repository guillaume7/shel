classdef WaterlevelInitialFactory
    methods(Static)
        function field = create(type, M, N, params)
            switch type
                case 'constant'
                    field = WaterlevelInitialFunctions.constantField(M, N, params.eta0);
                case 'gaussianBump'
                    field = WaterlevelInitialFunctions.gaussianBump(M, N, params.x0, params.y0, params.sx, params.sy, params.amplitude, params.dx, params.dy);
                case 'geostrophic'
                    field = WaterlevelInitialFunctions.geostrophicWaterLevel(M, N, params.eta0, params.gradient);
                case 'null'
                    field = WaterlevelInitialFunctions.nullField(M, N);
                otherwise
                    error('Unknown waterlevel initial condition type');
            end
        end
    end
end
