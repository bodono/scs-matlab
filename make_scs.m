gpu = false; % compile the gpu version of SCS
float = false; % using single precision (rather than double) floating points
int = false; % use 32 bit integers for indexing
% WARNING: OPENMP WITH MATLAB CAN CAUSE ERRORS AND CRASH, USE WITH CAUTION:
% openmp parallelizes the matrix multiply for the indirect solver (using CG):
flags.COMPILE_WITH_OPENMP = false;

flags.BLASLIB = '-lmwblas -lmwlapack';
% MATLAB_MEX_FILE env variable sets blasint to ptrdiff_t
flags.LCFLAG = '-DMATLAB_MEX_FILE -DUSE_LAPACK -DCTRLC=1 -DCOPYAMATRIX';
flags.INCS = '';
flags.LOCS = '';

common_scs = 'scs/src/linalg.c scs/src/cones.c scs/src/accel.c scs/src/cs.c scs/src/util.c scs/src/scs.c scs/src/ctrlc.c scs/src/normalize.c scs/src/scs_version.c scs/linsys/common.c scs_mex.c';
if (~isempty (strfind (computer, '64')))
    flags.arr = '-largeArrayDims';
else
    flags.arr = '';
end

if ( isunix && ~ismac )
    flags.link = '-lm -lut -lrt';
elseif  ( ismac )
    flags.link = '-lm -lut';
else
    flags.link = '-lut';
    flags.LCFLAG = sprintf('-DNOBLASSUFFIX %s', flags.LCFLAG);
end

if (float)
    flags.LCFLAG = sprintf('-DSFLOAT %s', flags.LCFLAG);
end
if (int)
    flags.INT = '';
else
    flags.INT = '-DDLONG';
end

if (flags.COMPILE_WITH_OPENMP)
    flags.link = strcat(flags.link, ' -lgomp');
end

compile_direct(flags, common_scs);
compile_indirect(flags, common_scs);
if (gpu)
    compile_gpu(flags, common_scs);
end

% compile scs_version
mex -O -Iscs/include scs/src/scs_version.c scs_version_mex.c -output scs_version

%%
clear data cones
disp('Example run:');
randn('seed',9)
m = 9;
n = 3;
data.A = sparse(randn(m,n));
data.b = randn(m,1);
data.c = randn(n,1);
cones.l = m;
[x,y,s,info] = scs_indirect(data,cones,[]);
[x,y,s,info] = scs_direct(data,cones,[]);


if (gpu)
    [x,y,s,info] = scs_gpu(data,cones,[]);
end

% test-warm start with solution
disp('Warm-starting:')
data.x = x;
data.y = y;
data.s = s;
[x,y,s,info] = scs_indirect(data,cones,[]);


disp('SUCCESSFULLY INSTALLED SCS')
disp('(If using SCS with CVX, note that SCS only supports CVX v3.0 or later).')
