classdef Grid
    methods (Static)
        function state = setup(state)
            % Setup grid coordinates and masks
            state.x0 = state.dx/2;
            state.y0 = state.dy/2;
            state.x = (1:state.M)' * (1:state.N);
            state.y = (1:state.M)' * (1:state.N);
            % Add more grid setup as needed
        end
    end
end
