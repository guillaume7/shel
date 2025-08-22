classdef MomentumForcing
    methods(Static)
        function drag = bottom_drag(u, H, Cd)
            % Compute bottom drag term for momentum equation
            % u: velocity (u or v)
            % H: water column height
            % Cd: bottom drag coefficient (scalar or matrix)
            % drag: bottom drag term (same size as u)
            drag = Cd .* u ./ max(H, 1e-6); % avoid division by zero
        end
    end
end
