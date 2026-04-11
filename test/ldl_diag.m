classdef ldl_diag < matlab.unittest.TestCase
    % Test the LDL extraction and C-level solve algorithm in pure MATLAB.

    methods (Test)
        function test_ldl_solve_symmetric(testCase)
            % Verify the C-level LDL extraction and solve algorithm.
            % Uses a symmetrized KKT matrix (matching our C code).
            rng(1234)
            n = 3; m = 9;
            A = sparse(randn(m, n));
            R_x = 0.1 * ones(n, 1);
            R_y = 0.1 * ones(m, 1);

            % Full symmetric KKT (what our C code builds after symmetrization)
            K = sparse([diag(R_x), A'; A, -diag(R_y)]);

            [L, D, perm] = ldl(K, 'vector');

            % Verify: K(perm,perm) = L*D*L'
            K_perm = K(perm, perm);
            testCase.verifyEqual(L*D*L', K_perm, 'AbsTol', 1e-12)

            % Replicate C-level solve (extract factors, forward/diag/backward)
            b = randn(n + m, 1);

            L_nodiag = L - speye(n + m);
            D_diag = full(diag(D, 0));
            D_sub = full(diag(D, -1));

            bp = b(perm);

            % Forward solve: (L_nodiag + I) * y = bp
            for i = 1:(n+m)
                for j = find(L_nodiag(:, i))'
                    bp(j) = bp(j) - L_nodiag(j, i) * bp(i);
                end
            end

            % Block diagonal solve
            i = 1;
            while i <= n + m
                if i < n + m && D_sub(i) ~= 0
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

            % Backward solve: (L_nodiag + I)' * x = z
            for i = (n+m):-1:1
                val = bp(i);
                for j = find(L_nodiag(:, i))'
                    val = val - L_nodiag(j, i) * bp(j);
                end
                bp(i) = val;
            end

            % Inverse permute
            x_ldl = zeros(n+m, 1);
            x_ldl(perm) = bp;

            % Reference: backslash
            x_ref = K \ b;

            testCase.verifyEqual(x_ldl, x_ref, 'AbsTol', 1e-10)
        end

        function test_ldl_requires_symmetric(testCase)
            % Verify that ldl on upper-tri-only gives WRONG results.
            % This documents why our C code must symmetrize before calling ldl.
            rng(5678)
            n = 3; m = 9;
            A = sparse(randn(m, n));
            R_x = 0.1 * ones(n, 1);
            R_y = 0.1 * ones(m, 1);

            K_full = sparse([diag(R_x), A'; A, -diag(R_y)]);
            K_upper = triu(K_full);

            b = randn(n + m, 1);

            % ldl on full symmetric: correct
            [L1, D1, p1] = ldl(K_full, 'vector');
            x1 = zeros(n+m, 1);
            y1 = L1' \ (D1 \ (L1 \ b(p1)));
            x1(p1) = y1;

            % ldl on upper-triangular only: incorrect
            [L2, D2, p2] = ldl(K_upper, 'vector');
            x2 = zeros(n+m, 1);
            y2 = L2' \ (D2 \ (L2 \ b(p2)));
            x2(p2) = y2;

            x_ref = K_full \ b;

            % Full symmetric should match
            testCase.verifyEqual(x1, x_ref, 'AbsTol', 1e-10)

            % Upper-tri-only should NOT match (documents the limitation)
            testCase.verifyGreaterThan(norm(x2 - x_ref), 0.1)
        end
    end
end
