classdef MomentumInitialFactory
    methods(Static)
        function field = create(type, M, N, value)
            switch type
                case 'constant'
                    field = MomentumInitialFunctions.constantField(M, N, value);
                case 'null'
                    field = MomentumInitialFunctions.nullField(M, N);
                otherwise
                    error('Unknown momentum initial condition type');
            end
        end
    end
end
