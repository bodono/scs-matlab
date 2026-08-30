function cmd = gpu_mex_command(flags, common_scs)
%GPU_MEX_COMMAND Construct the optional GPU indirect MEX build command.

flags.link = sprintf('-lcudart -lcublas -lcusparse %s', flags.link);
flags.INCS = sprintf('-I/usr/local/cuda/include %s', flags.INCS);
if ismac()
    flags.link = sprintf('-L/usr/local/cuda/lib %s', flags.link);
else
    % TODO probably not right for windows
    flags.link = sprintf('-L/usr/local/cuda/lib64 %s', flags.link);
end

% Match the native scs Makefile's gpu_indirect target: scs.c needs the
% INDIRECT definition, while the backend needs both its private solver and
% the shared GPU matrix implementation. Use the same integer-width flag as
% every other MEX backend so scs_int and cuSPARSE indices remain consistent.
cmd = sprintf(['mex -O -v %s %s %s %s -DINDIRECT=1 ' ...
    'COMPFLAGS="$COMPFLAGS %s" CFLAGS="$CFLAGS %s" ' ...
    'scs/linsys/gpu/indirect/private.c scs/linsys/gpu/gpu.c %s ' ...
    '-Iscs -Iscs/linsys -Iscs/include %s %s %s ' ...
    '-output matlab/scs_gpu'], ...
    flags.arr, flags.LCFLAG, flags.INT, common_scs, flags.COMPFLAGS, ...
    flags.CFLAGS, flags.INCS, flags.link, flags.LOCS, flags.BLASLIB);
