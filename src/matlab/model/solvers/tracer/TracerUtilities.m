classdef TracerUtilities
    methods (Static)
        function RHSt = ComputeSpaceT_UP_U(H, u, mask_u, Tr, K_L, dx_L)
            % Upwind scheme for tracer advection-diffusion
            [M,N]=size(H);
            Hm = zeros(M+1,N);
            TrUm = zeros(M+1,N);
            TrUp = zeros(M+1,N);
            TrDif = zeros(M+1,N);
            TrUm(1:M,:) = .5 * (abs(u(1:M,:)) - u(1:M,:)) .* Tr;
            TrUp(2:M+1,:) = .5 * (abs(u(2:M+1,:)) + u(2:M+1,:)) .* Tr;
            TrDif(2:M,:) = Tr(1:M-1,:) - Tr(2:M,:);
            Hm(2:M,:) = .5 * ( H(1:M-1,:) + H(2:M,:) );
            TrUm(M+1,:) = .5 * (abs(u(M+1,:)) - u(M+1,:)) * 0.;
            TrUp(1,:) = .5 * (abs(u(1,:)) + u(1,:)) * 0.;
            TrDif(1,:) = - Tr(1,:);
            TrDif(M+1,:) = Tr(M,:);
            Hm(1:M:M+1,:) = H(1:M-1:M,:);
            FFluxU = TracerUtilities.ComputeFaceFluxT_UP(Hm, TrUp, TrUm, TrDif, dx_L, K_L);
            RHSt = mask_u(1:M,:) .* FFluxU(1:M,:) - mask_u(2:M+1,:) .* FFluxU(2:M+1,:);
        end
        function RHSt = ComputeFaceFluxT_UP(Hm, TrUp, TrUm, TrDif, dx_L, K_L)
            RHSt = Hm .* (TrUp - TrUm + K_L / dx_L * TrDif) / dx_L;
        end
    end
end
