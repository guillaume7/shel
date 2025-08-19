classdef BottomInitialFactory
    methods(Static)
        function field = create(type, M, N, params)
            switch type
                case 'constant'
                    field = InitialConditionFunctions.constantField(M, N, params.d0);
                case 'step'
                    field = BottomInitialFunctions.stepField(M, N, params.d0, params.d0_step, params.stepIndex);
                case 'island'
                    field = BottomInitialFunctions.islandField(M, N, params.x0, params.y0, params.r, params.landValue, params.dx, params.dy);
                otherwise
                    error('Unknown bottom initial condition type');
            end
        end
    end
end
