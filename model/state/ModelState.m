classdef ModelState
    properties
        % Physical parameters
        dx, dy, dt, M, N, g, f, rho0, K
        Cd % Bottom drag coefficient
        rho_air % Air density
        uwind % Wind velocity (u-component)
        vwind % Wind velocity (v-component)
        % Fields
        eta, eta_old, eta_new
        H, H_old, H_new
        u, u_old, u_new, u_a
        v, v_old, v_new, v_a
        Tr, Tr_old, Tr_new
        mask, mask_u, mask_v, masknan
        p_surf % Atmospheric pressure deviation field (Pa)
        % Diagnostics
        ke, pe, eke, curl_t, potvorticity_t
        strechrate_t, shearrate_t, divergence_t
        sqstrain_t, sqstrechrate_t, sqshearrate_t
        enstrophy_t, okuboweiss_t
        % Time-dependent global properties
        time, volume, iKe, iPe, vtime, iVort, iEnst, iSqStrech, iSqShear, iSqStrain, iOWeiss, iMomentumU, iMomentumV, iEke
    end
    methods
        function obj = ModelState(params)
            % Initialize from params struct
            fields = fieldnames(params);
            for k = 1:numel(fields)
                obj.(fields{k}) = params.(fields{k});
            end
        end
    end
end
