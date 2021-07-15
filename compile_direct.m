function compile_direct(flags, common_scs)

if (flags.COMPILE_WITH_OPENMP)
    cmd = sprintf('mex -v -O %s %s %s %s COMPFLAGS="/openmp \\$COMPFLAGS" CFLAGS="\\$CFLAGS -fopenmp" -Iscs -Iscs/include %s', flags.arr, flags.LCFLAG, flags.INCS, flags.INT);
else
    cmd = sprintf ('mex -v -O %s %s %s %s -Iscs -Iscs/include -Iscs/linsys %s', flags.arr, flags.LCFLAG, flags.INCS, flags.INT);
end

amd_files = {'amd_order', 'amd_dump', 'amd_postorder', 'amd_post_tree', ...
    'amd_aat', 'amd_2', 'amd_1', 'amd_defaults', 'amd_control', ...
    'amd_info', 'amd_valid', 'amd_global', 'amd_preprocess', ...
    'SuiteSparse_config'} ;
for i = 1 : length (amd_files)
    cmd = sprintf ('%s scs/linsys/external/amd/%s.c', cmd, amd_files {i}) ;
end

cmd = sprintf ('%s scs/linsys/external/qdldl/qdldl.c %s scs/linsys/cpu/direct/private.c %s %s %s -output scs_direct', cmd, common_scs, flags.link, flags.LOCS, flags.BLASLIB);
eval(cmd);
