function compile_cudss(flags, common_scs)
%COMPILE_CUDSS Build the cuDSS GPU direct solver MEX (scs_cudss).
%
%   Requires the CUDA toolkit and the NVIDIA cuDSS library. Paths are
%   taken from the SCS_CUDA_PATH and SCS_CUDSS_PATH environment
%   variables, defaulting to /usr/local/cuda and /usr/local/cudss.

args = cudss_mex_command(flags, common_scs);
disp(strjoin(args, ' '));
mex(args{:});
