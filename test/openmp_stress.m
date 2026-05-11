classdef openmp_stress < matlab.unittest.TestCase
    % Stress test for the OpenMP-parallelized code paths in SCS.
    %
    % SCS has ``#pragma omp parallel for`` regions in
    % ``scs/linsys/scs_matrix.c`` (sparse matvec inside the indirect CG
    % solver) and ``scs/src/cones.c`` (cone projections). When the MEX
    % is built with ``use_open_mp = true`` these regions actually run
    % across threads, and a race condition would manifest as either
    % non-deterministic results across repeated solves of the same
    % problem or as outright crashes / NaNs.
    %
    % This test is intentionally run in both serial and OpenMP CI
    % builds: in serial it's a redundant correctness check (the loops
    % are sequential), in OpenMP mode it exercises the libiomp5-backed
    % parallel regions. The regression-pin job is the OpenMP CI matrix
    % entry — if a future change re-introduces a dual-runtime conflict
    % (or a race condition) this test catches it.

    properties
        data
        cones
    end

    properties (TestParameter)
        solver = {'default', 'indirect'}
    end

    methods(TestMethodSetup)
        function setup_problem(testCase)
            % Moderate-sized mixed-cone LP/SOC problem. Multiple SOC
            % blocks exercise the parallel cone-projection loop, and
            % the indirect solver exercises the parallel sparse matvec
            % during its CG iterations.
            rng(2026)
            n = 80;
            nl = 50;                            % LP cone
            nq_blocks = [40, 40, 40];           % three SOCs
            nq = sum(nq_blocks);
            m = nl + nq;

            A = sparse(randn(m, n));
            c = randn(n, 1);

            % Construct b so the problem is provably primal-feasible:
            % pick any x, place s strictly in the interior of each
            % cone, then set b = A*x + s.
            x_feas = randn(n, 1);
            s_feas = zeros(m, 1);
            s_feas(1:nl) = 1.0;                 % LP cone interior
            idx = nl + 1;
            for sz = nq_blocks
                s_feas(idx) = 1.0;              % SOC: t > ||v|| with v = 0
                idx = idx + sz;
            end
            testCase.data = struct( ...
                'A', A, ...
                'b', A * x_feas + s_feas, ...
                'c', c);
            testCase.cones = struct('l', nl, 'q', nq_blocks);
        end
    end

    methods (Test)
        function test_repeated_solves_are_consistent(testCase, solver)
            % Solve the same problem many times in a row and verify
            % every run lands on the same solution. A race in the
            % parallel matvec or cone-projection regions shows up as
            % run-to-run divergence; a hard race shows up as NaN or as
            % a non-``solved`` status.
            pars = openmp_stress.solver_pars(solver);
            pars.verbose = 0;
            n_runs = 20;
            xs = zeros(numel(testCase.data.c), n_runs);
            ys = zeros(numel(testCase.data.b), n_runs);
            for k = 1:n_runs
                [x, y, ~, info] = scs(testCase.data, testCase.cones, pars);
                testCase.verifyEqual(info.status, 'solved', ...
                    sprintf('run %d did not solve', k));
                testCase.verifyTrue(all(isfinite(x)), ...
                    sprintf('run %d produced non-finite primal', k));
                testCase.verifyTrue(all(isfinite(y)), ...
                    sprintf('run %d produced non-finite dual', k));
                xs(:, k) = x;
                ys(:, k) = y;
            end
            % All runs must agree. The tolerance is loose enough to
            % absorb harmless ULP-level drift from non-associative
            % parallel reductions, but tight enough to catch any real
            % nondeterminism: a race that scrambles even a handful of
            % iterates produces drift on the order of the SCS solve
            % tolerance (1e-3) or worse.
            for k = 2:n_runs
                testCase.verifyLessThan( ...
                    norm(xs(:, k) - xs(:, 1), inf), 1e-6, ...
                    sprintf('primal at run %d drifted from run 1', k));
                testCase.verifyLessThan( ...
                    norm(ys(:, k) - ys(:, 1), inf), 1e-6, ...
                    sprintf('dual at run %d drifted from run 1', k));
            end
        end
    end

    methods (Static)
        function pars = solver_pars(solver)
            pars = struct();
            if strcmp(solver, 'indirect'), pars.use_indirect = true; end
        end
    end
end
