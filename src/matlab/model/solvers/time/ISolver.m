classdef (Abstract) ISolver
    methods (Abstract)
        state = step(obj, state);
    end
end
