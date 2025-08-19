classdef OutputManager
    methods (Static)
        function saveOutput(state, t)
            % Save output fields to disk (expand as needed)
            % Example: save eta, u, v, H, diagnostics
            filename = sprintf('output_step_%04d.mat', t);
            eta = state.eta; u = state.u; v = state.v; H = state.H;
            save(filename, 'eta', 'u', 'v', 'H');
        end
        function plotModel(state)
            % Visualization logic (expand as needed)
            imagesc(state.eta); colorbar; title('Water Level');
        end
        function printStats(state, t)
            % Print simulation statistics
            fprintf('Step %d: Volume = %g, KE = %g, PE = %g\n', t, state.volume(t), state.iKe(t), state.iPe(t));
        end
    end
end
