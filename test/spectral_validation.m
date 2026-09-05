classdef spectral_validation < matlab.unittest.TestCase

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
        function test_rejects_shorter_nuc_n(testCase)
            cones = testCase.cones;
            cones.nuc_m = [2, 2];
            cones.nuc_n = 2;

            testCase.verifyError(@() testCase.solve(cones), 'scs:invalidCone')
        end

        function test_rejects_longer_nuc_n(testCase)
            cones = testCase.cones;
            cones.nuc_m = 2;
            cones.nuc_n = [2, 2];

            testCase.verifyError(@() testCase.solve(cones), 'scs:invalidCone')
        end

        function test_rejects_one_sided_nuclear_cone(testCase)
            cones = testCase.cones;
            cones.nuc_m = 2;

            testCase.verifyError(@() testCase.solve(cones), 'scs:invalidCone')
        end

        function test_rejects_shorter_sl_k(testCase)
            cones = testCase.cones;
            cones.sl_n = [2, 2];
            cones.sl_k = 1;

            testCase.verifyError(@() testCase.solve(cones), 'scs:invalidCone')
        end

        function test_rejects_longer_sl_k(testCase)
            cones = testCase.cones;
            cones.sl_n = 2;
            cones.sl_k = [1, 1];

            testCase.verifyError(@() testCase.solve(cones), 'scs:invalidCone')
        end

        function test_rejects_one_sided_sum_largest_cone(testCase)
            cones = testCase.cones;
            cones.sl_k = 1;

            testCase.verifyError(@() testCase.solve(cones), 'scs:invalidCone')
        end

        % A sparse double passes mxIsDouble but mxGetPr exposes only the
        % stored entries while numel reports the full length: reading the
        % full length walks past the array. Must be rejected, not read.
        function test_rejects_sparse_spectral_arrays(testCase)
            cones = testCase.cones;
            cones.nuc_m = sparse([2, 0]);
            cones.nuc_n = [2, 2];

            testCase.verifyError(@() testCase.solve(cones), 'scs:invalidCone')
        end

        function test_rejects_sparse_generic_cone_array(testCase)
            cones = testCase.cones;
            cones.q = sparse([2, 0]);

            testCase.verifyError(@() testCase.solve(cones), 'scs:invalidCone')
        end

        function test_rejects_complex_spectral_arrays(testCase)
            cones = testCase.cones;
            cones.sl_n = [2 + 1i, 2];
            cones.sl_k = [1, 1];

            testCase.verifyError(@() testCase.solve(cones), 'scs:invalidCone')
        end

        function test_rejects_nonfinite_spectral_arrays(testCase)
            cones = testCase.cones;
            cones.nuc_m = Inf;
            cones.nuc_n = 2;

            testCase.verifyError(@() testCase.solve(cones), 'scs:invalidCone')
        end

        function test_rejects_fractional_nuc_m(testCase)
            cones = testCase.cones;
            cones.nuc_m = 2.5;
            cones.nuc_n = 2;

            testCase.verifyError(@() testCase.solve(cones), 'scs:invalidCone')
        end

        function test_rejects_matrix_shaped_cone_array(testCase)
            cones = testCase.cones;
            cones.nuc_m = [2, 2; 2, 2];
            cones.nuc_n = [2, 2; 2, 2];

            testCase.verifyError(@() testCase.solve(cones), 'scs:invalidCone')
        end
    end

    methods
        function solve(testCase, cones)
            scs(testCase.data, cones, struct('verbose', 0));
        end
    end
end
