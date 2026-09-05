classdef input_validation < matlab.unittest.TestCase
% Type-safety of the MEX parsing layer: mxGetPr returns NULL for
% non-double arrays and exposes only stored entries of sparse arrays,
% so every scalar/array read must validate first (see #46).

    properties
        data
        cones
    end

    methods(TestMethodSetup)
        function setup_problem(testCase)
            testCase.data.A = sparse(-1);
            testCase.data.b = -1;
            testCase.data.c = 1;
            testCase.cones.l = 1;
        end
    end

    methods (Test)
        function test_logical_b_rejected(testCase)
            % mxGetPr returns NULL for non-double arrays; a logical b
            % previously reached the native copy path.
            data = testCase.data;
            data.b = true;
            testCase.verifyError( ...
                @() scs(data, testCase.cones, struct('verbose', 0)), ...
                'scs:invalidData')
        end

        function test_single_c_rejected(testCase)
            data = testCase.data;
            data.c = single(1);
            testCase.verifyError( ...
                @() scs(data, testCase.cones, struct('verbose', 0)), ...
                'scs:invalidData')
        end

        function test_complex_b_rejected(testCase)
            data = testCase.data;
            data.b = -1 + 2i;
            testCase.verifyError( ...
                @() scs(data, testCase.cones, struct('verbose', 0)), ...
                'scs:invalidData')
        end

        function test_complex_sparse_A_rejected(testCase)
            data = testCase.data;
            data.A = sparse(-1 + 1i);
            testCase.verifyError( ...
                @() scs(data, testCase.cones, struct('verbose', 0)), ...
                'scs:invalidData')
        end

        function test_dimension_mismatch_rejected(testCase)
            % A is 1x1 but numel(b) = 2. The scs.m wrapper's
            % scs_prepare_data catches this first with the same id...
            data = testCase.data;
            data.b = [-1; 0];
            testCase.verifyError( ...
                @() scs(data, testCase.cones, struct('verbose', 0)), ...
                'scs:invalidData')
        end

        function test_dimension_mismatch_rejected_at_mex(testCase)
            % ...and the MEX guard must hold on its own for direct
            % callers that bypass the wrapper: without it the CSC copy
            % reads past A's index arrays.
            data = testCase.data;
            data.b = [-1; 0];
            testCase.verifyError( ...
                @() scs_direct(data, testCase.cones, struct('verbose', 0)), ...
                'scs:invalidData')
        end

        function test_logical_warm_start_degrades_safely(testCase)
            % invalid warm-start input must fall back to zeros and still
            % solve, never reach the native copy
            data = testCase.data;
            data.x = true;
            data.y = true;
            data.s = true;
            pars = struct('verbose', 0, 'warm_start', 1);
            [x, ~, ~, info] = scs(data, testCase.cones, pars);
            testCase.verifyEqual(info.status, 'solved')
            testCase.verifyEqual(x, 1, 'AbsTol', 1e-3)
        end

        function test_logical_verbose_is_accepted(testCase)
            % struct('verbose', true) is idiomatic MATLAB; the old
            % *mxGetPr(tmp) read NULL for logicals.
            [x, ~, ~, info] = scs(testCase.data, testCase.cones, ...
                                  struct('verbose', false));
            testCase.verifyEqual(info.status, 'solved')
            testCase.verifyEqual(x, 1, 'AbsTol', 1e-3)
        end

        function test_char_setting_rejected(testCase)
            pars = struct('verbose', 0, 'alpha', 'fast');
            testCase.verifyError( ...
                @() scs(testCase.data, testCase.cones, pars), ...
                'scs:invalidSettings')
        end

        function test_sparse_setting_rejected(testCase)
            pars = struct('verbose', 0, 'alpha', sparse(1.5));
            testCase.verifyError( ...
                @() scs(testCase.data, testCase.cones, pars), ...
                'scs:invalidSettings')
        end

        function test_infinite_int_setting_rejected(testCase)
            pars = struct('verbose', 0, 'max_iters', Inf);
            testCase.verifyError( ...
                @() scs(testCase.data, testCase.cones, pars), ...
                'scs:invalidSettings')
        end

        function test_one_sided_box_rejected(testCase)
            data = testCase.data;
            data.A = sparse([1; 1]);
            data.b = [1; 0];
            cones = struct('bl', -1);
            testCase.verifyError( ...
                @() scs(data, cones, struct('verbose', 0)), ...
                'scs:invalidCone')
        end

        function test_fractional_int_setting_rejected(testCase)
            pars = struct('verbose', 0, 'max_iters', 2.5);
            testCase.verifyError( ...
                @() scs(testCase.data, testCase.cones, pars), ...
                'scs:invalidSettings')
        end

        function test_fractional_cone_array_rejected(testCase)
            cones.q = 2.5;
            testCase.verifyError( ...
                @() scs(testCase.data, cones, struct('verbose', 0)), ...
                'scs:invalidCone')
        end

        function test_sparse_scalar_cone_rejected(testCase)
            cones.l = sparse(1);
            testCase.verifyError( ...
                @() scs(testCase.data, cones, struct('verbose', 0)), ...
                'scs:invalidCone')
        end

        function test_complex_scalar_cone_rejected(testCase)
            cones.l = 1 + 2i;
            testCase.verifyError( ...
                @() scs(testCase.data, cones, struct('verbose', 0)), ...
                'scs:invalidCone')
        end

        function test_negative_scalar_cone_rejected(testCase)
            cones.l = -1;
            testCase.verifyError( ...
                @() scs(testCase.data, cones, struct('verbose', 0)), ...
                'scs:invalidCone')
        end

        function test_sparse_box_bound_rejected(testCase)
            data = testCase.data;
            data.A = sparse([1; 1]);
            data.b = [1; 0];
            cones = struct('bl', sparse(-1), 'bu', 1);
            testCase.verifyError( ...
                @() scs(data, cones, struct('verbose', 0)), ...
                'scs:invalidCone')
        end

        function test_infinite_box_bounds_are_accepted(testCase)
            % Box bounds legitimately use -Inf lower / +Inf upper; a
            % well-posed problem over the box [-1, Inf) must solve.
            % Rows: t = 1 (fixed), s1 = x, so x >= -1 with min x.
            data.A = sparse([0; -1]);
            data.b = [1; 0];
            data.c = 1;
            cones = struct('bl', -1, 'bu', Inf);
            [x, ~, ~, info] = scs(data, cones, struct('verbose', 0));
            testCase.verifyEqual(info.status, 'solved')
            testCase.verifyEqual(x, -1, 'AbsTol', 1e-3)
        end

    end
end
