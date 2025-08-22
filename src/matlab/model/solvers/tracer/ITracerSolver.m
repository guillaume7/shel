classdef (Abstract) ITracerSolver
    methods (Abstract)
        state = step(obj, state);
    end
end
