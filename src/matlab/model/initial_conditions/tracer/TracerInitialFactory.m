classdef TracerInitialFactory
    methods(Static)
        function field = create(type, M, N, params)
            switch type
                case 'constant'
                    field = TracerInitialFunctions.constantField(M, N, params.tracer0);
                case 'null'
                    field = TracerInitialFunctions.nullField(M, N);
                otherwise
                    error('Unknown tracer initial condition type');
            end
        end
    end
end
