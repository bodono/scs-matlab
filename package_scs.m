function package_scs()
    % PACKAGE_SCS  Programmatically create the SCS MATLAB Toolbox (.mltbx).
    % This script is used by GitHub Actions to automate distribution.

    try
        % 1. Define project root and metadata
        projectRoot = pwd;
        % SCS GUID - keep this constant to allow updates
        identifier = '9fe2fe4d-165e-4d2d-8511-91bcef1583b0';

        % 2. Create ToolboxOptions object (Requires R2023a or later)
        opts = matlab.addons.toolbox.ToolboxOptions(projectRoot, identifier);

        % 3. Configure Metadata
        opts.ToolboxName = 'SCS';
        opts.ToolboxVersion = '3.2.5';
        opts.AuthorName = 'Brendan O''Donoghue';
        opts.AuthorEmail = 'bodonoghue85@gmail.com';
        opts.Summary = 'Splitting Conic Solver for MATLAB';
        opts.Description = sprintf(['SCS (splitting conic solver) is a numerical optimization package ' ...
            'for solving large-scale convex cone problems, based on the paper ' ...
            '"Conic Optimization via Operator Splitting and Homogeneous Self-Dual Embedding" ' ...
            'by O''Donoghue, Chu, Parikh, and Boyd.']);

        % 4. Configure Output
        opts.OutputFile = fullfile(projectRoot, 'SCS.mltbx');

        % 5. Filter Toolbox Files
        % By default, ToolboxOptions includes everything in projectRoot.
        % We exclude Git, CI, and internal build files.
        allFiles = opts.ToolboxFiles;
        
        excludeList = { ...
            fullfile(projectRoot, '.git'), ...
            fullfile(projectRoot, '.github'), ...
            fullfile(projectRoot, '.gitmodules'), ...
            fullfile(projectRoot, '.bumpversion.cfg'), ...
            fullfile(projectRoot, 'package_scs.m'), ...
            fullfile(projectRoot, 'test') ...
        };
        
        isExcluded = false(size(allFiles));
        for i = 1:length(excludeList)
            isExcluded = isExcluded | startsWith(allFiles, excludeList{i});
        end
        
        opts.ToolboxFiles = allFiles(~isExcluded);

        % 6. Perform Packaging
        fprintf('Packaging SCS Toolbox version %s...\n', opts.ToolboxVersion);
        matlab.addons.toolbox.packageToolbox(opts);
        fprintf('Successfully created %s\n', opts.OutputFile);

    catch ME
        fprintf('Error during packaging: %s\n', ME.message);
        % Ensure we exit with error code for CI
        exit(1);
    end
end
