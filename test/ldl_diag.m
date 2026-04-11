classdef ldl_diag < matlab.unittest.TestCase
    % Diagnostic tests for the MATLAB LDL factorization and solve.
    %
    % These tests verify that our C-level LDL extraction and solve algorithm
    % (strip unit diagonal from L, extract D main/sub-diagonals, forward/
    % block-diagonal/backward solve with permutation) produces correct results.
    %
    % They also document a critical MATLAB API pitfall: sparse ldl requires
    % structurally symmetric input, even though the docs say "only the upper
    % triangle is referenced". See test_ldl_upper_tri_behavior.

    methods (Test)
        function test_ldl_solve_symmetric(testCase)
            % Verify our C-level solve algorithm on a small KKT system.
            rng(1234)
            n = 3; m = 9;
            A = sparse(randn(m, n));
            R_x = 0.1 * ones(n, 1);
            R_y = 0.1 * ones(m, 1);
            K = sparse([diag(R_x), A'; A, -diag(R_y)]);
            nm = n + m;

            [L, D, perm] = ldl(K, 'vector');
            testCase.verifyEqual(L*D*L', K(perm, perm), 'AbsTol', 1e-12)

            % Replicate C-level solve: extract L without diagonal, D diag/sub
            b = randn(nm, 1);
            L_nodiag = L - speye(nm);
            D_diag = full(diag(D, 0));
            D_sub = full(diag(D, -1));
            bp = b(perm);

            % Forward solve: (L + I) * y = bp
            for i = 1:nm
                for j = find(L_nodiag(:, i))'
                    bp(j) = bp(j) - L_nodiag(j, i) * bp(i);
                end
            end

            % Block diagonal solve: D * z = y
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

            % Backward solve: (L + I)' * x = z
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

        function test_ldl_solve_with_2x2_blocks(testCase)
            % Larger KKT system more likely to produce 2x2 Bunch-Kaufman
            % blocks in D. Verifies our block-diagonal solve handles them.
            rng(9876)
            n = 10; m = 30;
            A = sparse(randn(m, n));
            R_x = 0.01 * ones(n, 1);
            R_y = 0.01 * ones(m, 1);
            K = sparse([diag(R_x), A'; A, -diag(R_y)]);
            nm = n + m;

            [L, D, perm] = ldl(K, 'vector');
            D_sub = full(diag(D, -1));
            n_2x2_blocks = sum(D_sub ~= 0);

            % Replicate C-level solve
            b = randn(nm, 1);
            x_ldl = ldl_diag.c_level_solve(K, b);
            x_ref = K \ b;
            testCase.verifyEqual(x_ldl, x_ref, 'AbsTol', 1e-10, ...
                sprintf('Solve should be accurate (D has %d 2x2 blocks)', ...
                n_2x2_blocks))
        end

        function test_ldl_solve_with_P_matrix(testCase)
            % KKT with P (QP structure): K = [P+R_x, A'; A, -R_y]
            rng(5555)
            n = 8; m = 20;
            A = sparse(randn(m, n));
            P_half = randn(n, n);
            P = sparse(P_half * P_half');
            R_x = 0.1 * ones(n, 1);
            R_y = 0.1 * ones(m, 1);
            K = sparse([P + diag(R_x), A'; A, -diag(R_y)]);
            nm = n + m;

            b = randn(nm, 1);
            x_ldl = ldl_diag.c_level_solve(K, b);
            x_ref = K \ b;
            testCase.verifyEqual(x_ldl, x_ref, 'AbsTol', 1e-10, ...
                'KKT with P matrix should solve correctly')
        end

        function test_ldl_solve_multiple_rhs(testCase)
            % Solve with many right-hand sides to stress-test the algorithm.
            rng(7777)
            n = 5; m = 15;
            A = sparse(randn(m, n));
            R_x = 0.1 * ones(n, 1);
            R_y = 0.1 * ones(m, 1);
            K = sparse([diag(R_x), A'; A, -diag(R_y)]);
            nm = n + m;

            for trial = 1:20
                b = randn(nm, 1);
                x_ldl = ldl_diag.c_level_solve(K, b);
                x_ref = K \ b;
                testCase.verifyEqual(x_ldl, x_ref, 'AbsTol', 1e-10, ...
                    sprintf('RHS trial %d should solve correctly', trial))
            end
        end

        function test_ldl_upper_tri_behavior(testCase)
            % Documents that MATLAB's sparse ldl SILENTLY produces wrong
            % factors when given a structurally asymmetric (upper-tri only)
            % input. This is the reason scs_to_mxsparse_symmetric exists.
            %
            % MATLAB docs say ldl "only references the upper triangle", but
            % this means it reads only upper-tri *values* from a structurally
            % symmetric matrix. If the sparsity pattern itself is asymmetric
            % (i.e., entry (i,j) exists but (j,i) does not), ldl produces
            % garbage with no error or warning.
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

            % Upper-tri: verify that L*D*L' does NOT reconstruct K_full(p,p)
            [L2, D2, p2] = ldl(K_upper, 'vector');
            recon = L2 * D2 * L2';
            err_vs_full = norm(recon - K_full(p2,p2), 'fro');

            % The reconstruction error should be large — ldl produced garbage
            testCase.verifyGreaterThan(err_vs_full, 1e-6, ...
                'ldl(triu(K)) should NOT correctly factor the symmetric K')
        end
    end

    methods (Static)
        function x = c_level_solve(K, b)
            % Replicate the C-level LDL solve algorithm in MATLAB.
            % This mirrors forward_solve, diag_solve, backward_solve in
            % private/matlab_linsys/private.c
            nm = size(K, 1);
            [L, D, perm] = ldl(K, 'vector');
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

            x = zeros(nm, 1);
            x(perm) = bp;
        end
    end
end
