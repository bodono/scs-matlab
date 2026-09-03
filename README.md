scs-matlab
===

[![Build Status](https://github.com/bodono/scs-matlab/actions/workflows/build.yml/badge.svg)](https://github.com/bodono/scs-matlab/actions/workflows/build.yml)
[![Documentation](https://img.shields.io/badge/docs-online-brightgreen?logo=read-the-docs&style=flat)](https://www.cvxgrp.org/scs/)

Matlab interface for [SCS](https://github.com/cvxgrp/scs) 3.0.0 and higher.
The full documentation is available [here](https://www.cvxgrp.org/scs/).

## Installation

### 1. Download Toolbox (.mltbx) - Recommended
The easiest way to install SCS is to download the pre-compiled MATLAB Toolbox file:
1. Go to the [SCS Releases](https://github.com/bodono/scs-matlab/releases) page.
2. Download the latest `SCS.mltbx` file.
3. Open the file in MATLAB (or double-click it) to install.

This version comes with pre-compiled binaries for Windows, Linux, and both
Intel and Apple Silicon Macs, so you don't need a C compiler.

#### Linux binary compatibility

The pre-compiled Linux MEX files support x86-64 systems with glibc 2.31 or
newer. This is the minimum supported runtime, not an upper limit on the
user's glibc version. CI enforces this baseline by rejecting any binary that
requires a symbol newer than `GLIBC_2.31` before it is packaged.

The binaries are built with MATLAB R2022b on a pinned Ubuntu 22.04 runner.
The symbol check preserves compatibility with Ubuntu 20.04 while using
a [GitHub-hosted runner that is still supported](https://docs.github.com/en/actions/reference/runners/github-hosted-runners).

Systems older than glibc 2.31 can build SCS from source on the target machine.

#### Linux release compatibility checklist

When changing the minimum supported Linux distribution or MATLAB release,
or upgrading the build runner:

- Review the oldest supported runtime and its glibc version. The runner pin
  is a separate choice and can change while retaining the same ABI baseline.
- If dropping older-system support, update `SCS_LINUX_GLIBC_BASELINE` in
  `.github/workflows/build.yml`, the minimum above, and the release notes together.
- Keep the ABI check enabled; raise its baseline deliberately rather than
  removing it. Rebuild and load-test the release on the selected baseline
  and a current supported Linux system.

### 2. Build from Source
If you are on an unsupported platform or prefer to build from source:
1. Clone the repository recursively:
   ```bash
   git clone --recursive https://github.com/bodono/scs-matlab.git
   ```
2. In MATLAB, run the installer:
   ```matlab
   cd <path/to/scs-matlab>
   make_scs
   ```

## Usage

### One-shot solve

```matlab
data.P = sparse([3., -1.; -1., 2.]);
data.A = sparse([-1., 1.; 1., 0.; 0., 1.]);
data.b = [-1; 0.3; -0.5];
data.c = [-1.; -1.];
cone.z = 1;
cone.l = 2;

[x, y, s, info] = scs(data, cone, settings);
```

To warm-start, add fields `x`, `y`, `s` to the `data` struct from a
previous solve.

### Workspace reuse

When solving a sequence of problems where only `b` and/or `c` change (e.g.,
MPC, parameter sweeps), use the workspace API to avoid re-factorizing:

```matlab
work = scs_init(data, cone, settings);    % factorize once

[x, y, s, info] = scs_solve(work);       % solve

scs_update(work, b_new, []);             % update b ([] = unchanged)
[x, y, s, info] = scs_solve(work);       % re-solve (no re-factorization)

warm.x = x; warm.y = y; warm.s = s;
[x, y, s, info] = scs_solve(work, warm); % warm-started re-solve

scs_finish(work);                        % free workspace
```

### Solver backends

By default SCS uses MATLAB's built-in sparse LDL factorization (MA57 under
the hood). Alternatives:

```matlab
settings.use_qdldl = true;      % bundled QDLDL sparse direct solver
settings.use_indirect = true;    % conjugate gradient (iterative)
settings.dense = true;           % dense Cholesky (for dense A)
settings.gpu = true;             % cuDSS GPU direct solver (source build)
```

The GPU backend is the NVIDIA cuDSS direct solver (the GPU indirect
solver was removed in favor of it). It is optional, Linux-only in
practice, and not included in the precompiled toolbox. To build it from
source, install the CUDA toolkit and cuDSS (set `SCS_CUDA_PATH` /
`SCS_CUDSS_PATH` if not in `/usr/local/cuda` and `/usr/local/cudss`),
then:

```matlab
setenv('SCS_BUILD_GPU', 'true');
make_scs
```

### Cones

The `cone` struct fields correspond to the cone types. See the
[cone documentation](https://www.cvxgrp.org/scs/api/cones.html) for details.

| Field | Cone |
|-------|------|
| `z`   | Zero (equality constraints) |
| `l`   | Non-negative orthant |
| `bl`, `bu` | Box |
| `q`   | Second-order |
| `s`   | Semidefinite |
| `cs`  | Complex semidefinite |
| `ep`  | Primal exponential |
| `ed`  | Dual exponential |
| `p`   | Power |
