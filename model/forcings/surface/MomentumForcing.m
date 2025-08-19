classdef MomentumForcing
    methods(Static)
        function wind = wind_stress(u_or_v, rho0, rho_air, windvel, component)
            % Compute wind stress term for momentum equation
            % u_or_v: velocity field (u or v)
            % rho0: water density
            % rho_air: air density
            % windvel: wind velocity (u or v component)
            % component: 'u' or 'v'
            % wind: wind stress term (same size as u_or_v)
            Cw = 1.3e-3; % Typical wind drag coefficient
            wind = zeros(size(u_or_v));
            if strcmp(component, 'u')
                wind = (rho_air / rho0) * Cw * windvel * ones(size(u_or_v));
            elseif strcmp(component, 'v')
                wind = (rho_air / rho0) * Cw * windvel * ones(size(u_or_v));
            end
        end
    end
end
