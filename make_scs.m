gpu = false; % compile the gpu version of SCS
float = false; % using single precision (rather than double) floating points
int = false; % use 32 bit integers for indexing
% WARNING: OPENMP WITH MATLAB CAN CAUSE ERRORS AND CRASH, USE WITH CAUTION:
% openmp parallelizes the matrix multiply for the indirect solver (using CG)
% and some cone projections.
use_open_mp = false;

flags.BLASLIB = '-lmwblas -lmwlapack';
% MATLAB_MEX_FILE env variable sets blasint to ptrdiff_t
flags.LCFLAG = '-DMATLAB_MEX_FILE -DUSE_LAPACK -DCTRLC=1 -DCOPYAMATRIX -DGPU_TRANSPOSE_MAT -DVERBOSITY=0';
flags.INCS = '';
flags.LOCS = '';

common_scs = 'scs/src/linalg.c scs/src/cones.c scs/src/aa.c scs/src/util.c scs/src/scs.c scs/src/ctrlc.c scs/src/normalize.c scs/src/scs_version.c scs/linsys/scs_matrix.c scs/linsys/csparse.c scs/src/rw.c scs_mex.c';
if (contains(computer, '64'))
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

flags.CFLAGS = '';
flags.COMPFLAGS = '';

if (use_open_mp)
    flags.link = strcat(flags.link, ' -lgomp');
    flags.CFLAGS = sprintf('%s /openmp', flags.CFLAGS);
    flags.COMPFLAGS = sprintf('%s -fopenmp', flags.COMPFLAGS);
end

% add c99 to handle qldl comments
if (~ismac() && isunix())
    flags.CFLAGS = sprintf('%s -std=c99', flags.CFLAGS);
end

compile_direct(flags, common_scs);
compile_indirect(flags, common_scs);
if (gpu)
    compile_gpu(flags, common_scs);
end

% compile scs_version
mex -v -O -Iscs/include -Iscs/linsys scs/src/scs_version.c scs_version_mex.c -output scs_version

addpath '.'

disp('SUCCESSFULLY INSTALLED SCS')
disp('(If using SCS with CVX, note that SCS only supports CVX v3.0 or later).')
