classdef MomentumInitialFunctions
    methods(Static)
        function field = constantField(M, N, value)
            field = InitialConditionFunctions.constantField(M, N, value);
        end
        function field = nullField(M, N)
            field = InitialConditionFunctions.nullField(M, N);
        end
    end
end
