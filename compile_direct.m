function compile_direct(flags, common_scs)

cmd = sprintf('mex -O -v %s %s %s %s COMPFLAGS="$COMPFLAGS %s" CFLAGS="$CFLAGS %s" -Iscs -Iscs/linsys -Iscs/include', flags.arr, flags.LCFLAG, flags.INCS, flags.INT, flags.COMPFLAGS, flags.CFLAGS);

amd_files = {'amd_order', 'amd_dump', 'amd_postorder', 'amd_post_tree', ...
    'amd_aat', 'amd_2', 'amd_1', 'amd_defaults', 'amd_control', ...
    'amd_info', 'amd_valid', 'amd_global', 'amd_preprocess', ...
    'SuiteSparse_config'} ;
for i = 1 : length (amd_files)
    cmd = sprintf ('%s scs/linsys/external/amd/%s.c', cmd, amd_files {i}) ;
end

cmd = sprintf ('%s %s scs/linsys/external/qdldl/qdldl.c scs/linsys/cpu/direct/private.c %s %s %s -output scs_direct', cmd, common_scs, flags.link, flags.LOCS, flags.BLASLIB);
eval(cmd);
