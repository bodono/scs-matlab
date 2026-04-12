function package_scs()
    % PACKAGE_SCS  Programmatically create the SCS MATLAB Toolbox (.mltbx).
    % This script is used by GitHub Actions to automate distribution.

    % Determine project root from script location to be robust to caller's PWD
    projectRoot = fileparts(mfilename('fullpath'));
    oldDir = cd(projectRoot);
    cleanupObj = onCleanup(@() cd(oldDir));

    try
        % SCS GUID - keep this constant to allow updates
        identifier = '9fe2fe4d-165e-4d2d-8511-91bcef1583b0';

        % 1. Create ToolboxOptions object (Requires R2023a or later)
        opts = matlab.addons.toolbox.ToolboxOptions(projectRoot, identifier);

        % 2. Configure Metadata
        opts.ToolboxName = 'SCS';
        opts.ToolboxVersion = '3.2.5';
        opts.AuthorName = 'Brendan O''Donoghue';
        opts.AuthorEmail = 'bodonoghue85@gmail.com';
        opts.Summary = 'Splitting Conic Solver for MATLAB';
        opts.Description = sprintf(['SCS (splitting conic solver) is a numerical optimization package ' ...
            'for solving large-scale convex cone problems, based on the paper ' ...
            '"Conic Optimization via Operator Splitting and Homogeneous Self-Dual Embedding" ' ...
            'by O''Donoghue, Chu, Parikh, and Boyd.']);

        % 3. Configure Compatibility
        % Minimum release is tied to the release used to build the MEX binaries in CI
        opts.MinimumMatlabRelease = 'R2022b';
        
        % Disable MATLAB Online as this toolbox relies on native MEX binaries.
        % We also drop support for Intel Macs (Maci64) as we only build for
        % Apple Silicon (Macarm64) in CI.
        opts.SupportedPlatforms.Win64 = true;
        opts.SupportedPlatforms.Glnxa64 = true;
        opts.SupportedPlatforms.Maci64 = false;
        opts.SupportedPlatforms.Macarm64 = true;
        opts.SupportedPlatforms.MatlabOnline = false;

        % 4. Configure Output and Path
        opts.OutputFile = fullfile(projectRoot, 'SCS.mltbx');
        opts.ToolboxMatlabPath = {projectRoot};

        % 5. Filter Toolbox Files
        % We exclude Git metadata, CI files, and internal build artifacts from the bundle.
        allFiles = opts.ToolboxFiles;
        
        excludeList = { ...
            fullfile(projectRoot, '.git'), ...
            fullfile(projectRoot, '.github'), ...
            fullfile(projectRoot, '.gitmodules'), ...
            fullfile(projectRoot, '.bumpversion.cfg'), ...
            fullfile(projectRoot, 'package_scs.m'), ...
            fullfile(projectRoot, 'test'), ...
            fullfile(projectRoot, 'scs', '.git'), ...
            fullfile(projectRoot, 'scs', '.github'), ...
            fullfile(projectRoot, 'scs', 'test'), ...
            fullfile(projectRoot, 'scs', 'docs'), ...
            fullfile(projectRoot, 'scs', 'cmake'), ...
            fullfile(projectRoot, 'scs', 'Makefile'), ...
            fullfile(projectRoot, 'scs', 'scs.mk'), ...
            fullfile(projectRoot, 'scs', 'CMakeLists.txt'), ...
            fullfile(projectRoot, 'scs', 'include'), ...
            fullfile(projectRoot, 'scs', 'src'), ...
            fullfile(projectRoot, 'scs', 'linsys'), ...
            fullfile(projectRoot, 'private'), ...
            fullfile(projectRoot, 'scs', 'test') ...
        };
        
        isExcluded = false(size(allFiles));
        for i = 1:length(excludeList)
            isExcluded = isExcluded | startsWith(allFiles, excludeList{i});
        end
        
        opts.ToolboxFiles = allFiles(~isExcluded);

        % 6. Perform Packaging
        fprintf('Packaging SCS Toolbox version %s for %s+...\n', ...
            opts.ToolboxVersion, opts.MinimumMatlabRelease);
        matlab.addons.toolbox.packageToolbox(opts);
        fprintf('Successfully created %s\n', opts.OutputFile);

    catch ME
        fprintf('Error during packaging: %s\n', ME.message);
        % Ensure we exit with error code for CI
        exit(1);
    end
end
