classdef ldl_diag < matlab.unittest.TestCase
    % Diagnostic tests for MATLAB's ldl behavior with upper-triangular input.

    methods (Test)
        function test_ldl_solve_symmetric(testCase)
            % Verify our C-level solve algorithm is correct on symmetric input.
            rng(1234)
            n = 3; m = 9;
            A = sparse(randn(m, n));
            R_x = 0.1 * ones(n, 1);
            R_y = 0.1 * ones(m, 1);
            K = sparse([diag(R_x), A'; A, -diag(R_y)]);
            nm = n + m;

            [L, D, perm] = ldl(K, 'vector');
            testCase.verifyEqual(L*D*L', K(perm, perm), 'AbsTol', 1e-12)

            % Replicate C-level solve
            b = randn(nm, 1);
            L_nodiag = L - speye(nm);
            D_diag = full(diag(D, 0));
            D_sub = full(diag(D, -1));
            bp = b(perm);

            % Forward solve
            for i = 1:nm
                for j = find(L_nodiag(:, i))'
                    bp(j) = bp(j) - L_nodiag(j, i) * bp(i);
                end
            end

            % Block diagonal solve
            i = 1;
            while i <= nm
                if i < nm && D_sub(i) ~= 0
                    a = D_diag(i); bv = D_sub(i); d = D_diag(i+1);
                    det_val = a*d - bv*bv;
                    y0 = bp(i); y1 = bp(i+1);
                    bp(i) = (d*y0 - bv*y1) / det_val;
                    bp(i+1) = (a*y1 - bv*y0) / det_val;
                    i = i + 2;
                else
                    bp(i) = bp(i) / D_diag(i);
                    i = i + 1;
                end
            end

            % Backward solve
            for i = nm:-1:1
                val = bp(i);
                for j = find(L_nodiag(:, i))'
                    val = val - L_nodiag(j, i) * bp(j);
                end
                bp(i) = val;
            end

            x_ldl = zeros(nm, 1);
            x_ldl(perm) = bp;
            x_ref = K \ b;
            testCase.verifyEqual(x_ldl, x_ref, 'AbsTol', 1e-10)
        end

        function test_ldl_upper_tri_behavior(testCase)
            % Investigate what ldl actually does with upper-triangular input.
            rng(5678)
            n = 3; m = 9;
            A = sparse(randn(m, n));
            R_x = 0.1 * ones(n, 1);
            R_y = 0.1 * ones(m, 1);
            nm = n + m;

            K_full = sparse([diag(R_x), A'; A, -diag(R_y)]);
            K_upper = triu(K_full);

            b = randn(nm, 1);

            % Full symmetric: should be correct
            [L1, D1, p1] = ldl(K_full, 'vector');
            x1 = zeros(nm, 1);
            y1 = L1' \ (D1 \ (L1 \ b(p1)));
            x1(p1) = y1;
            x_ref = K_full \ b;
            testCase.verifyEqual(x1, x_ref, 'AbsTol', 1e-10, ...
                'ldl on full symmetric K should match backslash')

            % Upper-tri: check what L*D*L' actually reconstructs
            [L2, D2, p2] = ldl(K_upper, 'vector');
            recon = L2 * D2 * L2';

            % Does L*D*L' equal K_upper(p,p) (factoring the asymmetric matrix)?
            err_vs_upper = norm(recon - K_upper(p2,p2), 'fro');
            % Or does it equal K_full(p,p) (correctly using upper tri)?
            err_vs_full = norm(recon - K_full(p2,p2), 'fro');

            fprintf('ldl(K_upper): ||L*D*L'' - K_upper(p,p)||_F = %e\n', err_vs_upper);
            fprintf('ldl(K_upper): ||L*D*L'' - K_full(p,p)||_F  = %e\n', err_vs_full);

            % Check if ldl issued a warning
            [~, warnId] = lastwarn;
            fprintf('Last warning ID: %s\n', warnId);

            % If ldl is factoring K_upper (not K_full), the solve is wrong
            x2 = zeros(nm, 1);
            y2 = L2' \ (D2 \ (L2 \ b(p2)));
            x2(p2) = y2;
            err_solve = norm(x2 - x_ref) / norm(x_ref);
            fprintf('ldl(K_upper) solve relative error: %e\n', err_solve);

            % Record what happened for visibility, pass if symmetric works
            testCase.verifyLessThan(norm(x1 - x_ref), 1e-10, ...
                'Full symmetric solve must be accurate')
        end
    end
end
