function compile_indirect(flags, common_scs)

% compile indirect
if (flags.COMPILE_WITH_OPENMP)
    cmd = sprintf('mex -v -O %s %s %s %s -DINDIRECT COMPFLAGS="/openmp \\$COMPFLAGS" CFLAGS="\\$CFLAGS -fopenmp" scs/linsys/cpu/indirect/private.c %s -Iscs -Iscs/include %s %s %s -output scs_indirect',  flags.arr, flags.LCFLAG, common_scs, flags.INCS, flags.link, flags.LOCS, flags.BLASLIB, flags.INT);
else
    cmd = sprintf('mex -v -O %s %s %s %s -DINDIRECT scs/linsys/cpu/indirect/private.c %s -Iscs -Iscs/include -Iscs/linsys %s %s %s -output scs_indirect',  flags.arr, flags.LCFLAG, common_scs, flags.INCS, flags.link, flags.LOCS, flags.BLASLIB, flags.INT);
end
eval(cmd);
