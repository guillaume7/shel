classdef MomentumUtilities
    methods (Static)
        function av = fouraverage_u(v, M, N)
            av =  .25 * (v(2:M,1:N) + v(2:M,2:N+1) + v(1:M-1,1:N) + v(1:M-1,2:N+1));
        end
        function av = twoaverage_u(H, M, N)
            av = .5 * (H(2:M,2:N-1) + H(1:M-1,2:N-1));
        end
    end
end
