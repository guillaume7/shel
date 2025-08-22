classdef InitialConditionFunctions
    methods(Static)
        function field = constantField(M, N, value)
            % Create a constant field of size MxN with given value
            field = value * ones(M, N);
        end
        function field = nullField(M, N)
            % Zero field
            field = zeros(M, N);
        end
    end
end
