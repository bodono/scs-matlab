function compile_gpu(flags, common_scs)

flags.link = sprintf('-lcudart -lcublas -lcusparse %s', flags.link);
flags.INCS = sprintf('-I/usr/local/cuda/include %s', flags.INCS);
if (ismac())
    flags.link = sprintf('-L/usr/local/cuda/lib %s', flags.link);
else
    % TODO probably not right for windows
    flags.link = sprintf('-L/usr/local/cuda/lib64 %s', flags.link);
end

% compile gpu
cmd = sprintf('mex -O %s %s %s COMPFLAGS="$COMPFLAGS %s" CFLAGS="$CFLAGS %s" scs/linsys/gpu/indirect/private.c %s -Iscs -Iscs/linsys -Iscs/include %s %s %s -output scs_gpu',  flags.arr, flags.LCFLAG, common_scs, flags.COMPFLAGS, flags.CFLAGS, flags.INCS, flags.link, flags.LOCS, flags.BLASLIB);
eval(cmd);
