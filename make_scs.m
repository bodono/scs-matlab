scs_root = fileparts(mfilename('fullpath'));
old_dir = cd(scs_root);
src_dir = fullfile(scs_root, 'src');
matlab_dir = fullfile(scs_root, 'matlab');

% Add subfolders to path so we can find compilation and MATLAB scripts
addpath(src_dir);
addpath(matlab_dir);

gpu = false; % compile the gpu version of SCS
float = false; % using single precision (rather than double) floating points
int = false; % use 32 bit integers for indexing
% WARNING: OPENMP WITH MATLAB CAN CAUSE ERRORS AND CRASH, USE WITH CAUTION:
% openmp parallelizes the matrix multiply for the indirect solver (using CG)
% and some cone projections.
use_open_mp = false;

flags.BLASLIB = '-lmwblas -lmwlapack';
% MATLAB_MEX_FILE env variable sets blasint to ptrdiff_t
flags.LCFLAG = '-DMATLAB_MEX_FILE -DUSE_LAPACK -DCTRLC=1 -DCOPYAMATRIX -DGPU_TRANSPOSE_MAT -DVERBOSITY=0 -DUSE_SPECTRAL_CONES';
flags.INCS = '';
flags.LOCS = '';

common_scs = [ ...
    'scs/src/linalg.c scs/src/cones.c scs/src/exp_cone.c scs/src/aa.c ' ...
    'scs/src/util.c scs/src/scs.c scs/src/ctrlc.c scs/src/normalize.c ' ...
    'scs/src/scs_version.c scs/linsys/scs_matrix.c scs/linsys/csparse.c ' ...
    'scs/src/rw.c ' ...
    'scs/src/spectral_cones/logdeterminant/log_cone_Newton.c ' ...
    'scs/src/spectral_cones/logdeterminant/log_cone_IPM.c ' ...
    'scs/src/spectral_cones/logdeterminant/log_cone_wrapper.c ' ...
    'scs/src/spectral_cones/logdeterminant/logdet_cone.c ' ...
    'scs/src/spectral_cones/nuclear/ell1_cone.c ' ...
    'scs/src/spectral_cones/nuclear/nuclear_cone.c ' ...
    'scs/src/spectral_cones/sum-largest/sum_largest_cone.c ' ...
    'scs/src/spectral_cones/sum-largest/sum_largest_eval_cone.c ' ...
    'scs/src/spectral_cones/util_spectral_cones.c ' ...
    'src/scs_mex.c'];

if contains(computer, '64')
    flags.arr = '-largeArrayDims';
else
    flags.arr = '';
end

if isunix && ~ismac
    flags.link = '-lm -lut -lrt';
elseif ismac
    flags.link = '-lm -lut';
else
    flags.link = '-lut';
    flags.LCFLAG = ['-DNOBLASSUFFIX ' flags.LCFLAG];
end

if float
    flags.LCFLAG = ['-DSFLOAT ' flags.LCFLAG];
end

if int
    flags.INT = '';
else
    flags.INT = '-DDLONG';
end

flags.CFLAGS = '';
flags.COMPFLAGS = '';

if use_open_mp
    flags.link = [flags.link ' -lgomp'];
    flags.CFLAGS = [flags.CFLAGS ' /openmp'];
    flags.COMPFLAGS = [flags.COMPFLAGS ' -fopenmp'];
end

% add c99 to handle qldl comments
if isunix && ~ismac
    flags.CFLAGS = [flags.CFLAGS ' -std=c99'];
end

compile_direct(flags, common_scs);
compile_indirect(flags, common_scs);
compile_dense(flags, common_scs);
compile_matlab_direct(flags, common_scs);

if gpu
    compile_gpu(flags, common_scs);
end

% compile scs_version
cmd = sprintf('mex -v -O -I%s -I%s %s %s -output matlab/scs_version', ...
    fullfile('scs', 'include'), fullfile('scs', 'linsys'), ...
    fullfile('scs', 'src', 'scs_version.c'), ...
    fullfile('src', 'scs_version_mex.c'));
eval(cmd);

% Keep internal build helpers off the user's permanent MATLAB path.
rmpath(src_dir);

% Add to path and save
addpath(scs_root);
addpath(matlab_dir);
if savepath ~= 0
    fprintf('\n');
    fprintf('NOTE: Could not save the MATLAB path permanently.\n');
    fprintf('To use SCS in future sessions, add these lines to your startup.m:\n');
    fprintf('  addpath(''%s'');\n', scs_root);
    fprintf('  addpath(''%s'');\n', matlab_dir);
    fprintf('\n');
end

cd(old_dir);

disp('SUCCESSFULLY INSTALLED SCS')
disp('(If using SCS with CVX, note that SCS only supports CVX v3.0 or later).')
