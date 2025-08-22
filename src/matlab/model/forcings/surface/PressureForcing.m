classdef PressureForcing
    methods(Static)
        function dH = pressure_deviation(p_surf, rho0, g)
            % Compute water level adjustment due to atmospheric pressure deviation
            % p_surf: atmospheric pressure deviation field (Pa)
            % rho0: water density (kg/m^3)
            % g: gravity (m/s^2)
            % dH: water level adjustment (m)
            dH = p_surf ./ (rho0 * g);
        end
    end
end
