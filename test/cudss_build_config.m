classdef cudss_build_config < matlab.unittest.TestCase

    methods (Test)
        function test_cudss_command_shape(testCase)
            repo_root = fileparts(fileparts(mfilename('fullpath')));
            src_dir = fullfile(repo_root, 'src');
            addpath(src_dir);
            cleanup = onCleanup(@() rmpath(src_dir)); %#ok<NASGU>

            flags = struct( ...
                'arr', 'ARR_FLAG', ...
                'LCFLAG', 'LANGUAGE_FLAGS', ...
                'INT', 'INTEGER_WIDTH_FLAG', ...
                'COMPFLAGS', 'COMPILER_FLAGS', ...
                'CFLAGS', 'C_FLAGS', ...
                'INCS', 'INCLUDE_FLAGS', ...
                'link', 'LINK_FLAGS', ...
                'LOCS', 'LIBRARY_PATH_FLAGS', ...
                'BLASLIB', 'BLAS_FLAGS');

            args = cudss_mex_command(flags, 'COMMON_SCS_SOURCES');
            testCase.verifyTrue(iscellstr(args)); %#ok<ISCLSTR>
            cmd = strjoin(args, ' ');

            testCase.verifyEqual( ...
                count(cmd, 'scs/linsys/cudss/direct/private.c'), 1);
            testCase.verifyTrue(contains(cmd, '-lcudss'));
            testCase.verifyTrue(contains(cmd, '-lcudart'));
            % cuDSS supports only 32-bit SCS integers: the global
            % integer-width flag must be dropped from this backend's
            % command (a DLONG default must not leak through).
            testCase.verifyFalse(contains(cmd, 'INTEGER_WIDTH_FLAG'));
            testCase.verifyFalse(contains(cmd, '-DDLONG'));
            testCase.verifyTrue(contains(cmd, 'COMMON_SCS_SOURCES'));
            testCase.verifyTrue(contains(cmd, 'INCLUDE_FLAGS'));
            % fullfile uses the platform separator, so assert on the
            % argument cells rather than the joined string
            out_idx = find(strcmp(args, '-output'));
            testCase.verifyNumElements(out_idx, 1);
            testCase.verifyEqual(args{out_idx + 1}, ...
                                 fullfile('matlab', 'scs_cudss'));
            % direct backend: scs.c must NOT be built with INDIRECT
            testCase.verifyFalse(contains(cmd, '-DINDIRECT'));
            % Linux: rpath baked in, and cuBLAS made a direct dependency
            % (libcudss needs it but carries no RUNPATH of its own; GNU
            % DT_RUNPATH is not inherited by child resolution). Windows
            % resolves DLLs via PATH; macOS has no CUDA.
            expect_linux = isunix && ~ismac;
            testCase.verifyEqual(contains(cmd, '-Wl,-rpath'), expect_linux);
            testCase.verifyEqual(contains(cmd, '-lcublas'), expect_linux);
            testCase.verifyEqual(contains(cmd, '--no-as-needed'), expect_linux);
        end

        function test_gpu_pars_route_to_cudss(testCase)
            testCase.verifyEqual(scs_backend(struct('gpu', true)), 'scs_cudss');
            testCase.verifyEqual(scs_backend(struct('use_indirect', true)), 'scs_indirect');
            testCase.verifyEqual(scs_backend(struct()), 'scs_matlab_direct');
        end
    end
end
