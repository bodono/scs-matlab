function compile_gpu(flags, common_scs)
% compile gpu

cmd = gpu_mex_command(flags, common_scs);
disp(cmd);
eval(cmd);
