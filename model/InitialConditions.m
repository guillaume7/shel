classdef InitialConditions
    methods (Static)
        function state = setup(state, params)
            % Setup initial conditions for all test cases and fields
            M = state.M; N = state.N; dx = state.dx; dy = state.dy;
            % Coordinates for T, U, V, W cells
            state.x0 = dx/2; state.y0 = dy/2;
            % T-cell
            state.x = zeros(M,N); state.y = zeros(M,N);
            for j = 1:N
                state.x(:,j) = state.x0 + ((1:M)'/j - 1) * dx;
            end
            for i = 1:M
                state.y(i,:) = state.y0 + ((1:N)/i - 1) * dy;
            end
            % U-cell
            state.x_u = zeros(M+1,N); state.y_u = zeros(M+1,N);
            for j = 1:N
                state.x_u(:,j) = state.x0 + ((1:M+1)'/j - 1) * dx - dx/2;
            end
            for i = 1:M
                state.y_u(i,:) = state.y0 + ((1:N)/i - 1) * dy;
            end
            % V-cell
            state.x_v = zeros(M,N+1); state.y_v = zeros(M,N+1);
            for j = 1:N
                state.x_v(:,j) = state.x0 + ((1:M)'/j - 1) * dx;
            end
            for i = 1:M
                state.y_v(i,:) = state.y0 + ((1:N+1)/i - 1) * dy - dy/2;
            end
            % W-cell
            state.x_w = zeros(M,N); state.y_w = zeros(M,N);
            for j = 1:N
                state.x_w(:,j) = state.x0 + ((1:M)'/j - 1) * dx - dx/2;
            end
            for i = 1:M
                state.y_w(i,:) = state.y0 + ((1:N)/i - 1) * dy - dy/2;
            end
            % Bathymetry and mask
            state.d = params.d0 * ones(M,N);
            if params.tc_taylor
                % Submarine mount logic (expand as needed)
            elseif params.loadbathymetry
                % Load bathymetry from file (expand as needed)
            end
            if params.step
                state.d = Utilities.makestep(state.d, params.d0_step);
            end
            if params.tc_isla
                state.d = Utilities.makeisland(state.d, M*dx*0.5-params.isla_x0, N*dy*0.5-params.isla_y0, params.isla_R);
            end
            % Masks
            state.mask = ones(M,N);
            state.mask_u = ones(M+1,N);
            state.mask_v = ones(M,N+1);
            % Initial water level and velocity fields
            state.eta_old = params.eta0 * ones(M,N);
            state.eta = state.eta_old;
            state.eta_new = state.eta_old;
            state.H_old = state.eta_old + state.d;
            state.H = state.H_old;
            state.H_new = state.H_old;
            % Initial velocities
            state.u_old = zeros(M+1,N);
            state.u = state.u_old;
            state.u_new = state.u_old;
            state.u_a = state.u_old;
            state.v_old = zeros(M,N+1);
            state.v = state.v_old;
            state.v_new = state.v_old;
            state.v_a = state.v_old;
            % Tracer
            state.Tr = zeros(M,N);
            % Z fields
            state.curl_w = zeros(M+1,N+1);
            state.potvorticity_w = zeros(M+1,N+1);
            state.enstrophy_w = zeros(M+1,N+1);
            state.shearrate_w = zeros(M+1,N+1);
            state.sqshearrate_w = zeros(M+1,N+1);
            state.okuboweiss_w = zeros(M+1,N+1);
            % Visualization masks
            state.masknan = ones(M,N);
            % Diagnostics
            state.ke = zeros(M,N); state.pe = zeros(M,N);
            state.curl_t = zeros(M,N); state.potvorticity_t = zeros(M,N);
            state.strechrate_t = zeros(M,N); state.shearrate_t = zeros(M,N);
            state.enstrophy_t = zeros(M,N); state.sqstrain_t = zeros(M,N);
            state.okuboweiss_t = zeros(M,N); state.sqshearrate_t = zeros(M,N);
            state.sqstrechrate_t = zeros(M,N); state.divergence_t = zeros(M,N);
            state.sqdivergence_t = zeros(M,N);
            state.gradux = zeros(M,N); state.graduy = zeros(M,N);
            state.gradvx = zeros(M,N); state.gradvy = zeros(M,N);
            state.eke = zeros(M,N);
            % One-dimensional diagnostics
            L = params.L;
            state.vtime = nan(1,L); state.iKe = nan(1,L); state.iPe = nan(1,L);
            state.iEke = nan(1,L); state.volume = nan(1,L); state.iVort = nan(1,L);
            state.iEnst = nan(1,L); state.iSqShear = nan(1,L); state.iSqStrech = nan(1,L);
            state.iSqStrain = nan(1,L); state.iOWeiss = nan(1,L); state.iMomentumU = nan(1,L); state.iMomentumV = nan(1,L);
        end
    end
end
