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
% OpenMP parallelizes the matrix multiply for the indirect solver (using CG)
% and some cone projections. When enabled, this script links the MEX
% against MATLAB's own libiomp5 (the OpenMP runtime MATLAB already loads
% for its internal parallelism) rather than the compiler's default
% OpenMP runtime — see the ``if use_open_mp`` block below for the
% rationale and the per-platform setup. Windows is not currently
% supported with this approach; the build raises a clear error on
% Windows if you flip this on.
use_open_mp = false;

% Allow non-interactive callers (CI, scripts) to override the build
% knobs by setting environment variables before invoking this file.
% Interactive users editing the file directly are unaffected: an unset
% env var leaves the value at its literal default above.
if strcmpi(getenv('SCS_BUILD_GPU'), 'true');    gpu = true;         end
if strcmpi(getenv('SCS_USE_FLOAT'), 'true');    float = true;       end
if strcmpi(getenv('SCS_USE_INT'), 'true');      int = true;         end
if strcmpi(getenv('SCS_USE_OPENMP'), 'true');   use_open_mp = true; end

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
    % MATLAB ships Intel's libiomp5 and already loads it into its
    % process for internal parallelism. Linking a second OpenMP runtime
    % into the MEX (libgomp from gcc/MinGW, vcomp from MSVC) is
    % documented as undefined behavior — common failure modes are hangs,
    % wrong numerics, and crashes on MEX teardown. The fix is to compile
    % the MEX with the platform compiler's OpenMP frontend (so
    % ``#pragma omp`` regions are processed) but link against MATLAB's
    % own libiomp5 at link time. One OpenMP runtime in the process.
    matlab_arch = lower(computer('arch'));
    matlab_bin = fullfile(matlabroot, 'bin', matlab_arch);
    if ispc
        % MinGW's ``-fopenmp`` emits libgomp ABI calls (``GOMP_*``);
        % MATLAB's bundled ``libiomp5md.dll`` on Windows does not export
        % the ``GOMP_*`` compatibility aliases that the Linux build of
        % libiomp5 does. MSVC's ``/openmp`` emits ``__vcomp*`` calls
        % which libiomp5md also doesn't satisfy. Until we ship our own
        % GOMP-aliased iomp5 (or change SCS to use a different threading
        % primitive on Windows), fail loud rather than build a MEX that
        % hangs or crashes inside MATLAB at runtime.
        error('scs:openmpUnsupportedOnWindows', ...
            ['OpenMP is not yet supported on Windows for the MATLAB ' ...
             'MEX build (see comments in make_scs.m). Set ' ...
             'use_open_mp = false to continue.']);
    elseif ismac
        % clang's OpenMP ABI is compatible with Intel's libiomp5 by
        % design (clang's libomp and libiomp5 share the LLVM/Intel
        % OpenMP runtime origin). ``-Xclang -fopenmp`` enables pragma
        % processing in the frontend without triggering the driver's
        % auto-link of libomp.
        flags.CFLAGS = [flags.CFLAGS ' -Xclang -fopenmp'];
    else
        % gcc: ``-fopenmp`` at compile so the pragmas are processed
        % and the generated code emits ``GOMP_*`` runtime calls. At
        % link time we deliberately do NOT pass ``-fopenmp`` (which
        % would auto-link libgomp) — we substitute MATLAB's libiomp5
        % (which exports ``GOMP_*`` aliases on Linux) so the
        % GCC-emitted calls resolve into it. One runtime, no conflict.
        flags.CFLAGS = [flags.CFLAGS ' -fopenmp'];
    end
    % Pass libiomp5 by absolute path rather than ``-L<dir> -liomp5``.
    % Two reasons we can't use the latter on this code path:
    %   (1) ``flags.link`` is appended raw into the ``eval(cmd)``
    %       string in compile_*.m. Tokens that contain commas (e.g.
    %       ``-Wl,...``) are interpreted by MATLAB's parser as
    %       function-style arg separators, which then fails on the
    %       unary-minus prefixed link flags ("Invalid use of
    %       operator"). That rules out ``-Wl,-rpath,...``.
    %   (2) mex's mexopts template extracts ``-l*`` tokens into a
    %       ``$LIBS`` placeholder and ``-L*`` tokens into a separate
    %       ``$LIBPATHS`` placeholder, and emits LIBS *before*
    %       LIBPATHS in the gcc invocation. ``-liomp5 ... -L<path>``
    %       in that order means GNU ld hits the ``-l`` before it has
    %       seen its search path and reports ``cannot find -liomp5``.
    % An absolute path passed as a bare arg sidesteps both: mex sees
    % a ``.so``/``.dylib`` file and forwards it to the linker as a
    % direct input, no search-order dependency. We probe a small
    % list of candidate filenames because MATLAB has historically
    % shipped either ``libiomp5.so`` or ``libiomp5.so.5`` (versioned
    % soname) depending on release.
    if ismac
        candidates = {'libiomp5.dylib'};
    else
        candidates = {'libiomp5.so', 'libiomp5.so.5'};
    end
    libiomp5_path = '';
    for k = 1:numel(candidates)
        candidate = fullfile(matlab_bin, candidates{k});
        if isfile(candidate)
            libiomp5_path = candidate;
            break;
        end
    end
    if isempty(libiomp5_path)
        error('scs:libiomp5NotFound', ...
            'libiomp5 not found under %s (tried: %s)', ...
            matlab_bin, strjoin(candidates, ', '));
    end
    flags.link = sprintf('%s %s', flags.link, libiomp5_path);
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
disp(cmd);
eval(cmd);

% Keep internal build helpers off the user's permanent MATLAB path.
rmpath(src_dir);

% Add to path and save
addpath(matlab_dir);
if savepath ~= 0
    fprintf('\n');
    fprintf('NOTE: Could not save the MATLAB path permanently.\n');
    fprintf('To use SCS in future sessions, add this line to your startup.m:\n');
    fprintf('  addpath(''%s'');\n', matlab_dir);
    fprintf('\n');
end

cd(old_dir);

disp('SUCCESSFULLY INSTALLED SCS')
disp('(If using SCS with CVX, note that SCS only supports CVX v3.0 or later).')
