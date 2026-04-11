classdef ldl_diag < matlab.unittest.TestCase
    % Diagnostic test: replicate the C-level LDL extraction and solve
    % in pure MATLAB to verify correctness of the algorithm.

    methods (Test)
        function test_ldl_solve_upper_tri(testCase)
            % Build a KKT-like upper-triangular sparse matrix
            rng(1234)
            n = 3; m = 9;
            A = sparse(randn(m, n));
            R_x = 0.1 * ones(n, 1);
            R_y = 0.1 * ones(m, 1);

            % Full symmetric KKT
            K_full = [diag(R_x), A'; A, -diag(R_y)];

            % Upper-triangular only (what our C code passes to ldl)
            K_upper = triu(K_full);

            % Call ldl on upper-triangular
            [L, D, perm] = ldl(K_upper, 'vector');

            % Verify: K_upper(perm,perm) = L*D*L'
            K_perm = K_full(perm, perm);
            testCase.verifyEqual(L*D*L', K_perm, 'AbsTol', 1e-12)

            % Now replicate C-level solve
            b = randn(n + m, 1);

            % Extract L without diagonal (C code strips unit diagonal)
            L_nodiag = L - speye(n + m);

            % Extract D diagonal and sub-diagonal
            D_diag = full(diag(D, 0));
            D_sub = full(diag(D, -1));

            % C-level solve: permute, forward, diag, backward, unpermute
            % 0-indexed perm in C, but we use 1-indexed here
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
                    % 2x2 block
                    a = D_diag(i); bv = D_sub(i); d = D_diag(i+1);
                    det = a*d - bv*bv;
                    y0 = bp(i); y1 = bp(i+1);
                    bp(i) = (d*y0 - bv*y1) / det;
                    bp(i+1) = (a*y1 - bv*y0) / det;
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
            x_ref = K_full \ b;

            testCase.verifyEqual(x_ldl, x_ref, 'AbsTol', 1e-10)
        end

        function test_ldl_vs_full_symmetric(testCase)
            % Check that ldl on upper-tri gives same result as on full
            rng(5678)
            n = 3; m = 9;
            A = sparse(randn(m, n));
            R_x = 0.1 * ones(n, 1);
            R_y = 0.1 * ones(m, 1);

            K_full = sparse([diag(R_x), A'; A, -diag(R_y)]);
            K_upper = triu(K_full);

            b = randn(n + m, 1);

            [L1, D1, p1] = ldl(K_full, 'vector');
            [L2, D2, p2] = ldl(K_upper, 'vector');

            x1 = zeros(n+m, 1);
            y1 = L1' \ (D1 \ (L1 \ b(p1)));
            x1(p1) = y1;

            x2 = zeros(n+m, 1);
            y2 = L2' \ (D2 \ (L2 \ b(p2)));
            x2(p2) = y2;

            x_ref = K_full \ b;

            testCase.verifyEqual(x1, x_ref, 'AbsTol', 1e-10)
            testCase.verifyEqual(x2, x_ref, 'AbsTol', 1e-10)
        end
    end
end
