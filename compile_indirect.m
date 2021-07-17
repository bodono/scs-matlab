function compile_indirect(flags, common_scs)
% compile indirect
cmd = sprintf('mex -O %s %s %s %s COMPFLAGS="$COMPFLAGS %s" CFLAGS="$CFLAGS %s" scs/linsys/cpu/indirect/private.c %s -Iscs -Iscs/linsys -Iscs/include %s %s %s -output scs_indirect', flags.arr, flags.LCFLAG, common_scs, flags.INCS, flags.COMPFLAGS, flags.CFLAGS, flags.link, flags.LOCS, flags.BLASLIB, flags.INT);
eval(cmd);
