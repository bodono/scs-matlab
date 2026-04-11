classdef benchmark < matlab.unittest.TestCase
    % Head-to-head benchmark: matlab_ldl (default) vs qdldl.
    % Runs both solvers on a range of problem sizes, sparsity patterns, and
    % problem types (LP and QP), then prints an aggregated summary table.

    methods (Test)
        function test_benchmark(testCase)
            ns = [10, 50, 100, 200];
            m_ratios = [3, 5];
            densities = [0.1, 0.5, 1.0];
            types = {'LP', 'QP'};

            n_configs = length(ns) * length(m_ratios) * length(densities) * length(types);
            results = struct( ...
                'type', cell(n_configs, 1), ...
                'n', 0, 'm', 0, 'density', 0, ...
                'iter_default', 0, 'iter_qdldl', 0, ...
                'setup_default', 0, 'setup_qdldl', 0, ...
                'solve_default', 0, 'solve_qdldl', 0);

            idx = 0;
            for ti = 1:length(types)
                for ni = 1:length(ns)
                    for mi = 1:length(m_ratios)
                        for di = 1:length(densities)
                            n = ns(ni);
                            m = m_ratios(mi) * n;
                            density = densities(di);
                            prob_type = types{ti};
                            idx = idx + 1;

                            % Deterministic seed per config
                            rng(n * 1000 + m * 10 + round(density * 100) + ti)

                            % Build problem with guaranteed feasibility
                            if density >= 1.0
                                A = sparse(randn(m, n));
                            else
                                A = sprandn(m, n, density);
                            end
                            data.A = A;
                            K.l = m;

                            % Primal feasible: b = A*x + s, s > 0
                            x_feas = randn(n, 1);
                            s_feas = ones(m, 1);
                            data.b = A * x_feas + s_feas;
                            % Dual feasible: c = A'*y, y > 0 (guarantees boundedness)
                            y_feas = ones(m, 1);
                            data.c = A' * y_feas;

                            if strcmp(prob_type, 'QP')
                                P_half = sprandn(n, n, min(density * 2, 1.0));
                                data.P = sparse(P_half * P_half' + 0.1 * speye(n));
                            end

                            % Solve with default (matlab_ldl)
                            pars = struct('verbose', 0);
                            [x1, ~, ~, info1] = scs(data, K, pars);

                            % Solve with qdldl
                            pars_q = struct('verbose', 0, 'use_qdldl', true);
                            [x2, ~, ~, info2] = scs(data, K, pars_q);

                            % Verify both solvers agree on status
                            testCase.verifyEqual(info1.status, info2.status, ...
                                sprintf('%s n=%d m=%d d=%.1f: status mismatch', ...
                                prob_type, n, m, density))
                            % Compare solutions when both solved
                            if strcmp(info1.status, 'solved') && strcmp(info2.status, 'solved')
                                testCase.verifyEqual(x1, x2, 'AbsTol', 1e-3, ...
                                    sprintf('%s n=%d m=%d d=%.1f: solutions differ', ...
                                    prob_type, n, m, density))
                            end

                            % Record results
                            results(idx).type = prob_type;
                            results(idx).n = n;
                            results(idx).m = m;
                            results(idx).density = density;
                            results(idx).iter_default = info1.iter;
                            results(idx).iter_qdldl = info2.iter;
                            results(idx).setup_default = info1.setup_time;
                            results(idx).setup_qdldl = info2.setup_time;
                            results(idx).solve_default = info1.solve_time;
                            results(idx).solve_qdldl = info2.solve_time;
                        end
                    end
                end
            end

            % Print per-problem table
            fprintf('\n')
            fprintf('%-4s %5s %5s %7s | %6s %6s | %10s %10s | %10s %10s\n', ...
                'Type', 'n', 'm', 'density', ...
                'it_def', 'it_qdl', ...
                'setup_def', 'setup_qdl', ...
                'solve_def', 'solve_qdl')
            fprintf('%s\n', repmat('-', 1, 85))
            for i = 1:n_configs
                r = results(i);
                fprintf('%-4s %5d %5d %7.1f | %6d %6d | %8.2f ms %8.2f ms | %8.2f ms %8.2f ms\n', ...
                    r.type, r.n, r.m, r.density, ...
                    r.iter_default, r.iter_qdldl, ...
                    r.setup_default, r.setup_qdldl, ...
                    r.solve_default, r.solve_qdldl)
            end

            % Aggregate summary
            iters_def = [results.iter_default];
            iters_qdl = [results.iter_qdldl];
            setup_def = [results.setup_default];
            setup_qdl = [results.setup_qdldl];
            solve_def = [results.solve_default];
            solve_qdl = [results.solve_qdldl];
            total_def = setup_def + solve_def;
            total_qdl = setup_qdl + solve_qdl;

            fprintf('\n')
            fprintf('%-20s %12s %12s\n', '', 'matlab_ldl', 'qdldl')
            fprintf('%s\n', repmat('-', 1, 46))
            fprintf('%-20s %10d %12d\n', 'Median iters', ...
                round(median(iters_def)), round(median(iters_qdl)))
            fprintf('%-20s %10d %12d\n', 'Max iters', ...
                max(iters_def), max(iters_qdl))
            fprintf('%-20s %8.2f ms %10.2f ms\n', 'Median setup time', ...
                median(setup_def), median(setup_qdl))
            fprintf('%-20s %8.2f ms %10.2f ms\n', 'Max setup time', ...
                max(setup_def), max(setup_qdl))
            fprintf('%-20s %8.2f ms %10.2f ms\n', 'Median solve time', ...
                median(solve_def), median(solve_qdl))
            fprintf('%-20s %8.2f ms %10.2f ms\n', 'Max solve time', ...
                max(solve_def), max(solve_qdl))
            fprintf('%-20s %8.2f ms %10.2f ms\n', 'Median total time', ...
                median(total_def), median(total_qdl))
            fprintf('%-20s %8.2f ms %10.2f ms\n', 'Max total time', ...
                max(total_def), max(total_qdl))
            fprintf('\n')
        end
    end
end
