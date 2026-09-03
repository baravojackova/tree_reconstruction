%% START ============================================================
%  TreeQSM - reconstruction of tree structure from a point cloud
%  All file names are defined in ONE place (USER SETTINGS below).
%  Run the script cell by cell (Ctrl+Enter in each %% section).
%  ============================================================
clear
clc
%  ------------------------------------------------------------
%  USER SETTINGS - directory
%  ------------------------------------------------------------
% Folder with the TreeQSM source code (contains treeqsm.m, +myfun, ...)
src_dir = 'C:\Users\Spravce\Documents\BARA\01_Skeny_Babice\tree_reconstruction\scripts\matlab';

% Folder with the point clouds and where all results will be written
data_dir = 'C:\Users\Spravce\Documents\BARA\01_Skeny_Babice\tree_reconstruction\scripts\matlab';
%
%% 1) CLEAN START (OPTIONAL) ------------------------------------
%     CLEAN START (OPTIONAL) - delete ALL previously generated outputs
%     WHAT this block does: wipes every file this script has ever produced
%     in data_dir, for EVERY tree_id and EVERY run_tag at once (not just
%     the ones set above in section 1) - useful when you want to start
%     from a completely empty output database instead of accumulating
%     old runs.
%     WHY it needs two separate switches and a typed confirmation: this
%     is a MASS-DELETE of an entire folder's worth of results, and it is
%     very easy to accidentally run a whole script "cell by cell" and
%     land on this cell without meaning to. The switches below make sure
%     that (a) simply reaching this cell does nothing unless you have
%     deliberately edited clean_start, and (b) even with clean_start
%     enabled, nothing is actually deleted until you also disable the
%     dry run AND type a literal confirmation word.
%  ------------------------------------------------------------

% Switch 1 (master switch): must be hand-edited to true, otherwise this
% whole cell is a no-op. This is the main safety net against accidentally
% running this cell - the default (false) guarantees nothing happens.
clean_start = true;

% Switch 2 (dry run): with clean_start = true, this decides whether files
% are only LISTED (true, the safe default) or actually MOVED/DELETED
% (false, depending on clean_delete_permanently below). Keep this true the
% first time you enable clean_start, so you can review exactly what would
% happen before committing to it.
clean_dry_run = false;

% Switch 3 (archive vs. delete): false (default, SAFE) = ARCHIVE - move
% matched files into archive\<timestamp>\ instead of deleting them, fully
% recoverable. true = actually DELETE permanently (the old behaviour) -
% use this for genuine trial-run cleanups you want gone for good; leave it
% false for a "starting a new tree" cleanup where you'd rather have a
% safety net.
clean_delete_permanently = false;

if ~clean_start
    % Master switch is off - do nothing at all, not even list files.
    % This keeps a normal "run every cell in order" pass completely inert.
    fprintf('Clean start is disabled (clean_start = false) - skipping.\n');
else
    % ---- 1. Build the list of file-name PATTERNS that this script can
    % ever produce, across ALL tree_id/run_tag combinations (that's why
    % these use '*' wildcards instead of the specific names built in
    % section 2 below, which only cover the CURRENT tree_id/run_tag).
    clean_patterns = { ...
        '*_res_*.mat', ...              % first AND second run results (both contain '_res_')
        'volumes_*.mat', ...
        'volumes_*.csv', ...
        'dbh_*.txt', ...
        'height_*.txt', ...
        'taper_*.txt', ...
        'trunklen_*.txt', ...
        'trunklen_filtered_*.txt', ...
        'branchlen_*.txt', ...
        'branchlen_filtered_*.txt', ...
        'params_*.csv', ...
        'geom_*.txt', ...
        '*.mat', ...                    % catches saved point clouds too, e.g. "IND07_083.mat"
        };

    % ---- 2. Scan for every pattern above and collect matches, skipping
    % duplicates (a file can match more than one pattern, e.g. '*_res_*.mat'
    % and the broad '*.mat').
    %
    % WHY two search folders: make_models()/make_models_parallel() (the
    % TreeQSM functions called in section 12/15 below) do NOT save the
    % '*_res_*.mat' result files directly into data_dir - they hard-code a
    % "results" SUBFOLDER relative to the current directory (see
    % make_models.m / make_models_parallel.m: str = ['results/', savename]).
    % Since this script does cd(data_dir) in section 2, that means the
    % res_*.mat files actually live in [data_dir '\results'], not in
    % data_dir itself. Without also scanning that subfolder, this clean-up
    % section would silently miss every result file and report success
    % while leaving them all behind.
    clean_search_dirs = {data_dir, fullfile(data_dir, 'results')};

    clean_files = struct('name', {}, 'folder', {}, 'bytes', {});
    for d = 1:numel(clean_search_dirs)
        search_dir = clean_search_dirs{d};
        if ~isfolder(search_dir)
            continue   % e.g. no "results" subfolder yet - nothing to scan there
        end
        for p = 1:numel(clean_patterns)
            pattern = clean_patterns{p};
            matches = dir(fullfile(search_dir, pattern));
            for m = 1:numel(matches)
                f = matches(m);
                if f.isdir
                    continue   % dir() wildcards should only ever hit files here, but skip folders just in case
                end

                % ---- 3. SAFETY GUARD: never delete anything with "Cloud" in
                % its name, regardless of which pattern matched it. This is
                % the protection for raw input point clouds - both the
                % literal cloud_txt input file (e.g. "IND07_083 - Cloud_
                % branchbase_03.txt") and any hypothetical .mat file that
                % might one day be named similarly. Case-insensitive so
                % "cloud", "Cloud", "CLOUD" etc. are all protected the same way.
                if ~isempty(regexpi(f.name, 'Cloud', 'once'))
                    continue
                end

                % Skip if this exact file (by name+folder) is already in the list.
                already_listed = false;
                for e = 1:numel(clean_files)
                    if strcmp(clean_files(e).name, f.name) && strcmp(clean_files(e).folder, f.folder)
                        already_listed = true;
                        break
                    end
                end
                if ~already_listed
                    % Assign field-by-field (not "clean_files(end+1) = f") because
                    % dir() structs carry more fields (date, datenum, isdir, ...)
                    % than the 3-field placeholder clean_files was initialized
                    % with above - a whole-struct assignment between structs with
                    % different fields errors out ("dissimilar structures").
                    idx = numel(clean_files) + 1;
                    clean_files(idx).name   = f.name;   %#ok<SAGROW>  % small list, growth cost is irrelevant here
                    clean_files(idx).folder = f.folder;
                    clean_files(idx).bytes  = f.bytes;
                end
            end
        end
    end

    % ---- 4. Always print the full list of what would be (or was) deleted,
    % in both dry-run and real-delete mode, so the user can see exactly
    % what happened either way.
    fprintf('Clean start: %d matching file(s) found in %s\n', numel(clean_files), data_dir);
    total_bytes = 0;
    for i = 1:numel(clean_files)
        fprintf('  %-50s %10.1f KB\n', clean_files(i).name, clean_files(i).bytes / 1024);
        total_bytes = total_bytes + clean_files(i).bytes;
    end
    fprintf('Total: %d file(s), %.2f MB\n', numel(clean_files), total_bytes / (1024*1024));

    % archive_dir: only actually needed (and only actually CREATED, via
    % mkdir, in the real-move branch below) when clean_delete_permanently
    % is false - computed here regardless so the dry-run preview text can
    % name the exact folder that WOULD be used, without creating it.
    archive_dir = fullfile(data_dir, 'archive', datestr(now, 'yyyy-mm-dd_HHMM'));

    if clean_dry_run
        % ---- 5a. Dry run - list only, never touch the disk. Preview text
        % describes whichever mode (archive/delete) is actually active, so
        % it can never be mistaken for the other mode's behaviour.
        if clean_delete_permanently
            fprintf('DRY RUN - nothing deleted. Set clean_dry_run = false to actually PERMANENTLY DELETE.\n');
        else
            fprintf('DRY RUN - nothing archived. Set clean_dry_run = false to actually MOVE these files to %s\n', archive_dir);
        end
    else
        % ---- 5b. Real archive/delete - require a typed confirmation
        % first, so a stray re-run of this cell (e.g. Ctrl+Enter twice)
        % cannot touch anything without the user explicitly typing the
        % word again. The word itself (and what it does) matches whichever
        % mode is active, so Bara can't mistake one for the other right
        % before confirming.
        if clean_delete_permanently
            confirmation = input('Type SMAZAT (all caps) to PERMANENTLY DELETE the files listed above: ', 's');
            expected = 'SMAZAT';
        else
            confirmation = input('Type ARCHIVOVAT (all caps) to MOVE the files listed above to the archive folder: ', 's');
            expected = 'ARCHIVOVAT';
        end
        if ~strcmp(confirmation, expected)
            fprintf('Confirmation not received - nothing was touched.\n');
        else
            if ~clean_delete_permanently
                % Real move (not dry run) - create the archive folder now,
                % right before the first file actually needs it.
                mkdir(archive_dir);
            end
            for i = 1:numel(clean_files)
                file_path = fullfile(clean_files(i).folder, clean_files(i).name);
                try
                    if clean_delete_permanently
                        delete(file_path);
                        fprintf('Deleted: %s\n', file_path);
                    else
                        movefile(file_path, fullfile(archive_dir, clean_files(i).name));
                        fprintf('Archived to %s: %s\n', archive_dir, clean_files(i).name);
                    end
                catch ME
                    % try/catch per file so one locked/missing file cannot
                    % abort the rest of the list.
                    if clean_delete_permanently
                        fprintf('FAILED to delete %s (%s)\n', file_path, ME.message);
                    else
                        fprintf('FAILED to archive %s (%s)\n', file_path, ME.message);
                    end
                end
            end
        end
    end
end
%% 1B) USER SETTINGS ------------------------------------------
%  all manual setting
%  ------------------------------------------------------------
%
% run_tag is now AUTO-GENERATED from manual_patchdiam/man_PD1/man_PD2Min/
% man_PD2Max/simp_MaxOrder/simp_SmallRadii/simp_ReplaceIterations below
% (see section 2, "DERIVED NAMES") - not set here by hand any more, so a
% forgotten manual run_tag update can never silently overwrite a previous
% run's results with different settings. See section 2 for the exact
% auto-generation formula.
% --- tree identification -------------------------------------
tree_id   = 'IND01_054';           % short name used for ALL output files
cloud_txt = 'IND01_054.txt';   % input point cloud (text file, 3 columns X Y Z)

% --- number of models ----------------------------------------
n_models_first = 5;    % models per parameter combination, first (coarse) run
n_models_opt   = 25;   % models with the optimal inputs, second run

% false (default, SAFE) = if a res_file/res_new_file already exists for
% this EXACT run_tag, SKIP make_models and just load the existing file
% instead - protects against silently discarding an expensive
% reconstruction by accidentally re-running with identical settings.
% true = always recompute and overwrite, even if a matching file already
% exists.
overwrite_existing_reconstruction = false;

% --- parallel computing --------------------------------------
use_parallel = true;   % true  = make_models_parallel
                       % false = make_models (single core, fallback)
n_workers = 0;         % 0 = derive automatically from the number of tasks

% --- parameter ranges for define_input -----------------------
% define_input(P, nPD1, nPD2Min, nPD2Max) = how many values are tested
% for each of the three patch-diameter parameters
% !!! if manual input nPD = 1
nPD1    = 2;
nPD2Min = 2;
nPD2Max = 2;

% --- MANUAL PatchDiam (see section 9) ------------------------
manual_patchdiam = false;   % false = keep everything from define_input

% PatchDiam1 (rought first cover) has to be t ≥ PatchDiam2Max (gentle cover)

man_PD1    = 0.08;   % PatchDiam1     - AdQSM paper, Indonesian site 0,08
man_PD2Min = 0.02;   % PatchDiam2Min
man_PD2Max = 0.07;   % PatchDiam2Max
man_MinCylRad = 0.0025; % MinCylRad default 0.0025

% --- which model is simplified and exported ------------------
use_optimal = true;    % true  = use the optimal model from select_optimum
                       % false = use res.QSMs(model_index)
model_index = 1;       % used only when use_optimal = false

plot_optimal = true;   % true = plot the optimal QSM before simplification

% --- simplification setting  -----------------
% MaxOrder           Maximum branching order, higher order branches removed
% SmallRadii         Minimum acceptable radius for a branch at its base
% ReplaceIterations  Number of iterations for replacing two concecutive
%                    cylinders inside one branch with one longer cylinder
% --- simplification settings ---------------------------------
simp_MaxOrder          = 8;
simp_SmallRadii        = 0.005;
simp_ReplaceIterations = 2;
simp_Plot              = 1;
simp_Disp              = 1;

% Export the "Filtered <10cm" group (branches/trunk restricted to
% diameter >= 10 cm, matching how the destructive field reference was
% physically measured) into VolumeTable and into the
% trunklen_filtered_/branchlen_filtered_ text files. Set to false to
% skip this export entirely (e.g. when you don't need a same-cutoff
% comparison against the destructive reference for this run) - the
% diagnostic cylinder/volume printout in the console still runs either
% way, only the CSV/text-file EXPORT is affected.
export_filtered_10cm = true;

% Master switch for the ANSYS geometry export (section 20) - default OFF,
% same "explicit opt-in" pattern as clean_start/export_filtered_10cm above,
% so a normal cell-by-cell run does NOT write geom_*.txt files on every
% single iteration (section 20 used to run unconditionally). Set to true
% only when you actually want this run's geometry exported.
export_to_ansys = false;

% Optional override for the ANSYS export's filename tag - if left blank
% (''), the export falls back to whichever run's model is actually being
% exported (run_tag for whatever's live in memory, or
% EXPORT_FROM_SAVED_RUN_TAG below when that's set), so the geometry files
% are clearly labeled with whichever reconstruction variant produced them.
% Set to a non-blank string to use that instead (e.g. for a one-off export
% you want named something more memorable). See section 20 for where this
% is actually applied (as ansys_tag, not used directly).
ansys_export_name = '';

% Blank (default) = export whatever simplification is CURRENTLY LIVE in
% memory (today's behaviour, via the ansys_source switch in section 20).
% Non-blank = instead load 'simplified_<tree_id>_<THIS VALUE>.mat' from
% disk (see section 16c's unconditional save) and export THAT, ignoring
% ansys_source and whatever's currently live in memory - lets you
% retroactively export any previously-saved simplification variant, even
% one no longer in the workspace (e.g. after re-running section 2 + 16 for
% a second variant, per the "Solution A" workflow).
EXPORT_FROM_SAVED_RUN_TAG = '';

%% 2) DERIVED NAMES -------------------------------------------
%  - built automatically from tree_id + run_tag
%     Normally you do NOT edit this block
%  ------------------------------------------------------------
% run_tag: AUTO-GENERATED from the actual settings above (manual_patchdiam/
% man_PD1/man_PD2Min/man_PD2Max/simp_MaxOrder/simp_SmallRadii/
% simp_ReplaceIterations), computed HERE (not up in section 1c, where those
% settings are defined) so it always reflects whatever those variables are
% actually set to for THIS run - a forgotten manual run_tag edit after
% changing a setting can no longer silently overwrite a previous run's
% results under the same name. No manual version prefix (v1/v2/...) - not
% wanted; if you re-run with IDENTICAL settings, you'll get the identical
% tag and legitimately overwrite that run's own results, which is correct.
%
% Auto mode's actual PatchDiam values are only known DURING the run (from
% the auto-search inside define_input/treeqsm itself) - they intentionally
% do NOT appear in this tag, which is built from settings known BEFORE the
% run starts. They're captured separately, after the fact, in
% params_<tree>_<run>.csv (section 19) - that file remains the source of
% truth for auto mode's actual PD values, not this tag.
%
% Formula lives in compute_run_tag() (a local function at the end of this
% file, same pattern as find_disconnected_islands) - NOT inlined here -
% so section 16's sanity check (which needs the exact same formula, to
% catch "changed a setting but forgot to re-run this section") can call
% the identical code instead of a second, hand-copied formula that could
% drift out of sync with this one.
run_tag = compute_run_tag(manual_patchdiam, man_PD1, man_PD2Min, ...
    man_PD2Max, simp_MaxOrder, simp_SmallRadii, simp_ReplaceIterations);
fprintf('Auto-generated run_tag: %s\n', run_tag);

mat_name     = tree_id;                        % .mat file with the point cloud
res_name     = [tree_id '_res_'     run_tag];  % results of the first run
res_new_name = [tree_id '_res_new_' run_tag];  % results of the second run

mat_file     = [mat_name     '.mat'];
res_file     = [res_name     '.mat'];
res_new_file = [res_new_name '.mat'];

vol_file      = ['volumes_' tree_id '_' run_tag];       % volume table
export_prefix = ['geom_'    tree_id '_' run_tag '_'];   % geometry exports

% Make TreeQSM functions visible, then work inside the data folder
addpath(genpath(src_dir));
cd(data_dir);

fprintf('Tree: %s   Run: %s\n', tree_id, run_tag);
fprintf('Input cloud  : %s\n', cloud_txt);
fprintf('Point cloud  : %s\n', mat_file);
fprintf('Results 1    : %s\n', res_file);
fprintf('Results 2    : %s\n', res_new_file);
fprintf('Volume table : %s.mat / .csv\n', vol_file);
fprintf('Geometry     : %s*.txt\n', export_prefix);

%% 3) IMPORT  ------------------------------------------------
%  the point cloud from the text file
%  ------------------------------------------------------------
if ~isfile(cloud_txt)
    error('Input file not found: %s', fullfile(data_dir, cloud_txt));
end
P = load(cloud_txt);
fprintf('Loaded %d points.\n', size(P,1));

%% 4) SHIFT ---------------------------------------------------
%  the coordinate system so that the cloud starts at 0
%  ------------------------------------------------------------
P(:,1) = P(:,1) - min(P(:,1));   % X
P(:,2) = P(:,2) - min(P(:,2));   % Y
P(:,3) = P(:,3) - min(P(:,3));   % Z

fprintf('Extent: X %.2f m, Y %.2f m, Z %.2f m\n', ...
    max(P(:,1)), max(P(:,2)), max(P(:,3)));

%% 5) PLOT the point cloud ------------------------------------
%  for visual check
%  ------------------------------------------------------------
% Explicit, high figure number (10) - NOT auto-numbered - so it can never
% collide with figures 1/2, which simplify_qsm.m (a shared TreeQSM library
% function, left untouched) hard-codes internally for the optimal/
% simplified model plots below. A plain, unnumbered figure() call here
% used to get auto-assigned Figure 1 (being the first figure created),
% which simplify_qsm's hard-coded figure(1) call would then silently
% overwrite (no hold on, so its plot3 call clears the axes) - see
% CHANGELOG/investigation notes for the full trace.
% MATLAB does not allow a numeric figure handle together with
% Name/Value pairs in one figure(...) call ("Numeric figure handles not
% supported with parameter-value pairs") - so the number and the Name are
% set in two separate statements instead.
figure(10);
set(gcf, 'Name', ['Point cloud - ' tree_id]);
plot3(P(:,1), P(:,2), P(:,3), '.k', 'MarkerSize', 1);
axis equal;                      % same scale on all axes
grid off;
xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]');
title(['Cleaned and positioned tree: ' tree_id], 'Interpreter', 'none');

%% 6) SAVE the point cloud -----------------------------------
%  as a .mat file
%  ------------------------------------------------------------
save(mat_file, 'P');
fprintf('Point cloud saved as %s\n', mat_file);

%% 7) RE-ENTRY POINT ------------------------------------------
%   - start here if the .mat file already exists
%  ------------------------------------------------------------
load(mat_file);                  % loads variable P

%%  8) DEFINE the input parameters automatically----------------
%    define_input first runs create_input, so ALL other fields
%     (nmin1, TaperCor, ParentCor, ...) are already filled in.
%     It then sets only PatchDiam1/2Min/2Max and BallRad1/2.
%  ------------------------------------------------------------
inputs = define_input(P, nPD1, nPD2Min, nPD2Max);

inputs.name = tree_id;           % name used inside TreeQSM output files
inputs.tree = 1;
inputs.plot = 0;                 % 0 = no plots during the batch run
inputs.disp = 1;                 % 1 = short text output

% Keep a copy of the automatic values for comparison in section 9
auto.PD1    = inputs.PatchDiam1;
auto.PD2Min = inputs.PatchDiam2Min;
auto.PD2Max = inputs.PatchDiam2Max;
auto.BR1    = inputs.BallRad1;
auto.BR2    = inputs.BallRad2;
auto.MinCylRad = inputs.MinCylRad;

fprintf('\n--- define_input (automatic) ---\n');
fprintf('PatchDiam1    = %.4f    BallRad1 = %.4f\n', auto.PD1,    auto.BR1);
fprintf('PatchDiam2Min = %.4f\n',                    auto.PD2Min);
fprintf('PatchDiam2Max = %.4f    BallRad2 = %.4f\n', auto.PD2Max, auto.BR2);
fprintf('Implied stem radius Rstem = %.3f m (PatchDiam1 * 3)\n', auto.PD1*3);   %estimate how big stem radius is estimated by QSM

%% 9) MANUAL PatchDiam - BallRad derived with the SAME formulas------
%       that define_input uses:
%       BallRad1 = max( PD1    + 1.50*Res , min(1.25*PD1,    PD1    + 0.025) )
%       BallRad2 = max( PD2Max + 1.25*Res , min(1.20*PD2Max, PD2Max + 0.025) )
%     Res = point resolution, recovered from the automatic values.
%  ------------------------------------------------------------
if manual_patchdiam

    % --- recover the point resolution Res --------------------
    % If the automatic BallRad1 came from the "PD1 + 1.5*Res" branch,
    % then Res = (BallRad1 - PatchDiam1) / 1.5.
    Res = (auto.BR1 - auto.PD1)/1.5;

    % If the other branch was active, the value above is an upper
    % bound only. Cap it so it cannot inflate the manual BallRad.
    Res_cap = (min(1.25*auto.PD1, auto.PD1 + 0.025) - auto.PD1) / 1.5;
    if ~isfinite(Res) || Res < 0
        Res = 0;
    elseif Res <= Res_cap
        Res = 0;      % the min() branch was active -> Res not recoverable
    end
    fprintf('\nPoint resolution Res = %.4f m\n', Res);

    % --- apply the manual patch diameters --------------------
    inputs.PatchDiam1    = man_PD1;
    inputs.PatchDiam2Min = man_PD2Min;
    inputs.PatchDiam2Max = man_PD2Max;
    inputs.MinCylRad = man_MinCylRad;

    % --- derive the ball radii, same formulas as define_input -
    inputs.BallRad1 = max([inputs.PatchDiam1 + 1.5*Res, ...
        min(1.25*inputs.PatchDiam1, inputs.PatchDiam1 + 0.025)]);

    inputs.BallRad2 = max([inputs.PatchDiam2Max + 1.25*Res, ...
        min(1.20*inputs.PatchDiam2Max, inputs.PatchDiam2Max + 0.025)]);

    % --- comparison printout ---------------------------------
    fprintf('\n--- PatchDiam / BallRad ---\n');
    fprintf('                  auto     manual\n');
    fprintf('PatchDiam1      %6.4f    %6.4f\n', auto.PD1,    inputs.PatchDiam1);
    fprintf('BallRad1        %6.4f    %6.4f\n', auto.BR1,    inputs.BallRad1);
    fprintf('PatchDiam2Min   %6.4f    %6.4f\n', auto.PD2Min, inputs.PatchDiam2Min);
    fprintf('PatchDiam2Max   %6.4f    %6.4f\n', auto.PD2Max, inputs.PatchDiam2Max);
    fprintf('BallRad2        %6.4f    %6.4f\n', auto.BR2,    inputs.BallRad2);
    fprintf('BallRad1/PD1    %6.2f    %6.2f\n', ...
        auto.BR1/auto.PD1, inputs.BallRad1/inputs.PatchDiam1);
     fprintf('MinCylRad        %6.4f    %6.4f\n', auto.MinCylRad,    inputs.MinCylRad);

    % --- sanity check ----------------------------------------
    if inputs.BallRad1 <= inputs.PatchDiam1 || ...
       inputs.BallRad2 <= inputs.PatchDiam2Max
        warning('BallRad should be larger than the corresponding PatchDiam.');
    end
    if inputs.PatchDiam2Min >= inputs.PatchDiam2Max
        warning('PatchDiam2Min should be smaller than PatchDiam2Max.');
    end

else
    disp('Using PatchDiam and BallRad from define_input.');
end

%% 9b) START the parallel pool (run this BEFORE step 10)------
%       Doing it here separates pool problems from QSM problems
%  ------------------------------------------------------------
if use_parallel
    n_tasks = nPD1 * nPD2Min * nPD2Max * n_models_first;

    if n_workers == 0
        n_req = min(n_tasks, 32);   % never more workers than tasks
    else
        n_req = n_workers;
    end

    pool = gcp('nocreate');
    if isempty(pool) || pool.NumWorkers < n_req
        delete(pool);
        pool = parpool('Processes', n_req);
    end
    fprintf('Tasks: %d, pool: %d workers.\n', n_tasks, pool.NumWorkers);
end

%% 10) FIRST RUN - models over the parameter grid------------
%  
%  ------------------------------------------------------------
tic     % start time monitoring
if ~overwrite_existing_reconstruction && isfile(fullfile('results', res_file))
    fprintf('res_file already exists (%s) - skipping first-run reconstruction.\n', res_file);
    fprintf('(set overwrite_existing_reconstruction = true to force a fresh recompute)\n');
else
    if use_parallel
        QSMs = make_models_parallel(mat_name, res_name, n_models_first, inputs);
    else
        QSMs = make_models(mat_name, res_name, n_models_first, inputs);
    end
    fprintf('First QSM run finished in %.1f min.\n', toc/60);   % toc - stop time moniroting
end

%% 11) LOAD the results of the first run------------------------
%  
%  ------------------------------------------------------------
res = load(res_file);
fprintf('Loaded %d models from %s\n', numel(res.QSMs), res_file);

%% 12) OPTIMISATION - select the best models--------------
%  (based on point-to-cylinder distances)
%  ------------------------------------------------------------
[TreeData_O, OptModels, OptInputs, OptQSM] = select_optimum(res.QSMs);

%% 13) SECOND RUN - more models with the optimal parameters----
%  ------------------------------------------------------------
tic
if ~overwrite_existing_reconstruction && isfile(fullfile('results', res_new_file))
    fprintf('res_new_file already exists (%s) - skipping second-run reconstruction.\n', res_new_file);
    fprintf('(set overwrite_existing_reconstruction = true to force a fresh recompute)\n');
else
    if use_parallel
        QSMs_new = make_models_parallel(mat_name, res_new_name, n_models_opt, OptInputs);
    else
        QSMs_new = make_models(mat_name, res_new_name, n_models_opt, OptInputs);
    end
    fprintf('Second QSM run finished in %.1f min.\n', toc/60);
end

%% 14) LOAD the results of the second run---------------------
%  ------------------------------------------------------------
res_new = load(res_new_file);
fprintf('Loaded %d models from %s\n', numel(res_new.QSMs), res_new_file);

%% 15) PRECISION - combine both runs to get better std estimates---
%  ------------------------------------------------------------
[TreeData_E, OptQSMs_E, OptQSM_E] = estimate_precision( ...
    res.QSMs, res_new.QSMs, TreeData_O, OptModels);

%% 16) SELECT the source model, then SIMPLIFY it--------------
%  ------------------------------------------------------------
% Indices of the optimal group (needed here and in 16b)
%   OptModels{1} = indices of all models of the winning combination
%   OptModels{2} = index of the single representative model (OptQSM)
%
% MaxOrder      Maximum branching order, higher order branches removed
% SmallRadii    Minimum acceptable radius for a branch at its base
% ReplaceIterations Number of iterations for replacing two concecutive
%                     cylinders inside one branch with one longer cylinder
%
% --- simplification settings ---------------------------------
% simp_MaxOrder          = 8;
% simp_SmallRadii        = 0.005;
% simp_ReplaceIterations = 2;
% simp_Plot              = 1;
% simp_Disp              = 1;
%-------------------------------------------------------------
clear QSM_simple_clean   % avoid picking up a stale result from an
                          % earlier simplification pass in this same
                          % MATLAB session (see Solution A workflow -
                          % this section is re-run multiple times
                          % without clearing the whole workspace)

% ---- SAFETY CHECK: run_tag must match CURRENT settings ------------
% Catches "changed simp_*/manual_patchdiam/man_PD* but forgot to
% re-run section 2" - every downstream export in this pass (16c's
% simplified_file, volume table, dbh/height/params, geom_*.txt) would
% otherwise silently be written under the WRONG (stale) run_tag.
expected_run_tag = compute_run_tag(manual_patchdiam, man_PD1, man_PD2Min, ...
    man_PD2Max, simp_MaxOrder, simp_SmallRadii, simp_ReplaceIterations);
if ~strcmp(run_tag, expected_run_tag)
    error(['run_tag (''%s'') does not match what your CURRENT settings ' ...
           'would produce (''%s''). You changed manual_patchdiam/man_PD*/' ...
           'simp_* since section 2 last ran - re-run section ' ...
           '"2) DERIVED NAMES" before continuing.'], run_tag, expected_run_tag);
end

if iscell(OptModels)
    idx     = double(OptModels{1}(:))';
    idx_rep = double(OptModels{2}(1));
else
    idx     = double(OptModels(:))';
    idx_rep = idx(1);
end
fprintf('Optimal group: %d models, indices %s\n', numel(idx), mat2str(idx));
fprintf('Representative model: index %d\n', idx_rep);

if use_optimal
    QSM_opt = OptQSM(1);                 % optimal model from select_optimum
    disp('Source model: optimal QSM.');
else
    if model_index > numel(res.QSMs)
        error('model_index = %d, but only %d models exist.', ...
            model_index, numel(res.QSMs));
    end
    QSM_opt = res.QSMs(model_index);     % manually chosen model
    fprintf('Source model: res.QSMs(%d).\n', model_index);
end

% QSM_opt is kept EXACTLY as produced by select_optimum/estimate_precision,
% with no island cleaning applied - it is only a raw, unfiltered REFERENCE
% value used later (section 17, 'Optimal'/'Optimal (single)' rows) so you
% can compare against the model you actually work with. The model you
% actually use going forward (islands detected/removed, ANSYS export,
% etc.) is always QSM_simple, built below and cleaned in sections
% 15b/15c further down (moved there so island cleaning runs on the
% SIMPLIFIED model, which is what matters for the real workflow).
%
% simp_Plot (below) makes simplify_qsm.m draw TWO figures of its own:
% Figure 1 = the optimal (pre-simplification) model, i.e. QSM_opt as
%            plotted before simplify_qsm modifies it.
% Figure 2 = the simplified (post-simplification) model, i.e. QSM_simple.
% Both figure numbers are HARD-CODED inside simplify_qsm.m itself (a
% shared TreeQSM library function) - do NOT change them here. The point-
% cloud plot in section 5 above was moved to an explicit Figure 10 for
% exactly this reason (it used to auto-number as Figure 1 and get
% silently overwritten by simplify_qsm's Figure 1 call).
QSM_simple = simplify_qsm(QSM_opt, simp_MaxOrder, ...
    simp_SmallRadii, simp_ReplaceIterations, simp_Plot, simp_Disp);

%% 16b) DETECT DISCONNECTED BRANCH ISLANDS ------
%  - REVIEW BEFORE REMOVING
% This section only DETECTS and PLOTS potential disconnected branch
% fragments (cylinder.parent == 0 for a non-root cylinder) - it does
% NOT remove anything. Look at the resulting figure: islands are
% drawn in red on top of the tree in gray. Only run the NEXT section
% (16c) if you decide, after reviewing the plot, that these really
% are broken/noise fragments that should be removed.
% MOVED here (after simplify_qsm) and now reads from QSM_simple(end)
% instead of QSM_opt: the user's actual working model is the simplified
% one (fewer cylinders, used for ANSYS), so island cleaning should run on
% THAT model, not on the raw pre-simplification one. QSM_simple(end) is
% the same "the simplified model" choice already used elsewhere in this
% script (e.g. V_simp = vols(QSM_simple(end)) in section 17 below).
%  ------------------------------------------------------------

% cylinder.parent(k) is the row index (into the SAME cylinder table) of
% cylinder k's parent cylinder. A well-formed QSM has exactly ONE cylinder
% with parent == 0: cylinder 1, the trunk base (the true root, with no
% parent at all). If any OTHER cylinder also has parent == 0, that means
% TreeQSM/AdQSM lost track of how it connects to the rest of the tree -
% it and everything built on top of it form a separate "island" floating
% disconnected from the main structure. These usually come from noisy or
% incomplete point-cloud data and typically represent little volume, but
% should be reviewed visually (not blindly deleted) before removal.
parent_arr = QSM_simple(end).cylinder.parent;
radius_arr = QSM_simple(end).cylinder.radius;
length_arr = QSM_simple(end).cylinder.length;
start_arr  = QSM_simple(end).cylinder.start;   % (n_cyl,3) xyz of each cylinder's base
total_tree_volume = sum(pi .* radius_arr.^2 .* length_arr);

% find_disconnected_islands() is the local function defined at the very
% end of this file (see "LOCAL FUNCTIONS" there) - it walks the
% parent/child tree and returns one list of cylinder indices per island.
island_groups = find_disconnected_islands(parent_arr);

fprintf('\n--- Disconnected branch islands found: %d ---\n', numel(island_groups));
fprintf('Total tree volume: %.4f m3\n', total_tree_volume);

% NOTE: uses "island_idx" (NOT "idx") on purpose - section 16 above
% already assigned the script variable "idx" to the winning-combination
% MODEL indices, and section 17 below reads that same "idx" again
% (V_opt = vols(res.QSMs(idx))). Since this is a plain script, every
% section shares ONE workspace, so reusing "idx" here would silently
% overwrite it with cylinder indices instead and break section 17.
all_island_indices = [];
for oi = 1:numel(island_groups)
    island_idx = island_groups{oi};
    vol = sum(pi .* radius_arr(island_idx).^2 .* length_arr(island_idx));
    pct = vol / total_tree_volume * 100;
    z_range = [min(start_arr(island_idx,3)), max(start_arr(island_idx,3))];
    fprintf('  Island %d: %d cylinders, %.5f m3 (%.2f %% of tree), height %.2f-%.2f m\n', ...
        oi, numel(island_idx), vol, pct, z_range(1), z_range(2));
    all_island_indices = [all_island_indices, island_idx];
end

% Report the TOTAL volume across all islands combined (not just per
% island above), so you can see at a glance how much volume removing
% every island would cost - guarded with isempty() so sum() over an
% empty all_island_indices (zero islands found) does not run/print at all.
if ~isempty(island_groups)
    total_island_vol = sum(pi .* radius_arr(all_island_indices).^2 .* ...
        length_arr(all_island_indices));
    total_island_pct = total_island_vol / total_tree_volume * 100;
    fprintf('Total island volume: %.5f m3 (%.2f %% of tree)\n', ...
        total_island_vol, total_island_pct);
end

if isempty(island_groups)
    fprintf('No disconnected islands found.\n');
else
    % --- plot: whole tree in light gray, islands in red on top ---
    figure('Name', ['Disconnected branch islands - ' tree_id]);
    plot3(start_arr(:,1), start_arr(:,2), start_arr(:,3), '.', ...
          'Color', [0.3 0.3 0.3], 'MarkerSize', 6);
    hold on
    plot3(start_arr(all_island_indices,1), start_arr(all_island_indices,2), ...
          start_arr(all_island_indices,3), '.', 'Color', 'r', 'MarkerSize', 18);
    xlabel('x [m]'); ylabel('y [m]'); zlabel('z [m]');
    title(sprintf('%d disconnected island(s), %.2f %% of tree volume', ...
          numel(island_groups), sum(pi.*radius_arr(all_island_indices).^2 .* ...
          length_arr(all_island_indices))/total_tree_volume*100));
    axis equal
    grid off
    view(3)
    hold off
end

%% 16c) REMOVE DISCONNECTED ISLANDS-----------------------
%   (run manually, only after reviewing the plot above)
% Creates QSM_simple_clean - a COPY of QSM_simple(end) with the island
% cylinders removed and every remaining cylinder's "parent" index
% remapped to the new (shifted) row numbers. QSM_simple itself is left
% untouched, so re-running earlier sections is unaffected. Renamed from
% the former QSM_opt_clean, since this now cleans the SIMPLIFIED model,
% not the raw optimal one.
%  ------------------------------------------------------------

if isempty(island_groups)
    fprintf('No islands to remove - QSM_simple_clean not created.\n');
else
    n_cyl_total = numel(parent_arr);
    keep_mask = true(n_cyl_total, 1);
    keep_mask(all_island_indices) = false;   % mark every island cylinder for removal

    % old_to_new maps an OLD row index to its NEW row index after the
    % island rows are deleted (rows shift up to fill the gaps). A kept
    % cylinder that used to be row 37 might become row 30, for example -
    % every "parent" reference has to follow that same shift or it would
    % end up pointing at the wrong (unrelated) cylinder.
    old_to_new = zeros(n_cyl_total, 1);
    old_to_new(keep_mask) = 1:sum(keep_mask);

    QSM_simple_clean = QSM_simple(end);   % start from the simplified model, not QSM_opt
    cyl_fields = fieldnames(QSM_simple_clean.cylinder);
    for f = 1:numel(cyl_fields)
        field_name = cyl_fields{f};
        field_val = QSM_simple_clean.cylinder.(field_name);
        % Only touch fields that have one ROW per cylinder (size along
        % dim 1 equal to n_cyl_total) - e.g. radius/length/parent/start.
        % Anything else (a scalar setting, etc.) is left untouched.
        if size(field_val, 1) == n_cyl_total
            QSM_simple_clean.cylinder.(field_name) = field_val(keep_mask, :);
        end
    end

    % Remap parent indices: any parent > 0 (i.e. not itself a root) must
    % now point at the NEW row number of that same parent cylinder.
    new_parent = QSM_simple_clean.cylinder.parent;
    nonzero = new_parent > 0;
    new_parent(nonzero) = old_to_new(new_parent(nonzero));
    QSM_simple_clean.cylinder.parent = new_parent;

    r_c = QSM_simple_clean.cylinder.radius;
    L_c = QSM_simple_clean.cylinder.length;
    order_c = QSM_simple_clean.cylinder.BranchOrder;
    V_cyl_c = pi .* r_c.^2 .* L_c;

    total_vol_clean  = sum(V_cyl_c);
    stem_vol_clean   = sum(V_cyl_c(order_c == 0));
    branch_vol_clean = sum(V_cyl_c(order_c >= 1));
    V_rep_clean = [total_vol_clean, stem_vol_clean, branch_vol_clean];

    fprintf('QSM_simple_clean created: %d cylinders (was %d), removed %.5f m3 (%.2f %%).\n', ...
        sum(keep_mask), n_cyl_total, total_tree_volume - total_vol_clean, ...
        (total_tree_volume - total_vol_clean)/total_tree_volume*100);
end

% ---- PERSIST the final simplified model to disk, UNCONDITIONALLY -------
% Runs in BOTH branches above (islands found -> QSM_simple_clean, no
% islands -> QSM_simple_clean never created) and regardless of
% export_to_ansys (section 20) - every simplification pass gets saved
% automatically, so a second simplification variant (re-running section 2
% + jumping back to section 16 with different simp_* settings) never
% silently loses the FIRST variant's result from memory with no way to
% recover it. Section 20's EXPORT_FROM_SAVED_RUN_TAG can later reload any
% of these files by the run_tag that produced them.
if exist('QSM_simple_clean', 'var')
    QSM_final = QSM_simple_clean(1);
else
    QSM_final = QSM_simple(1);   % no islands - simplified alone IS the final result
end
simplified_file = ['simplified_' tree_id '_' run_tag '.mat'];
save(simplified_file, 'QSM_final');
fprintf('Simplified model saved to %s\n', simplified_file);

%% 16d) PLOT the simplified-------------------------------------
%   (no islands) model for visual comparison
%       against the raw "Simplified" model and the islands plot (15b).
%       Only meaningful when islands were actually found and removed
%       (guard mirrors the exist('QSM_simple_clean','var') check used
%       for the "Simplified (no islands)" row in the volume table, 17).
%  ------------------------------------------------------------

if isempty(island_groups)
    fprintf('No islands removed - skipping Simplified (no islands) plot.\n');
elseif ~exist('QSM_simple_clean', 'var')
    fprintf('QSM_simple_clean not found - skipping Simplified (no islands) plot.\n');
else
    fig_clean = figure('Name', ['Simplified model, no islands - ' tree_id]);
    if exist('plot_cylinder_model', 'file')
        % Built-in TreeQSM plotter - renders actual cylinders (not just points)
        plot_cylinder_model(QSM_simple_clean.cylinder, 'order', fig_clean.Number);
    else
        % Fallback: lightweight point plot, same style as the 15b islands plot
        start_clean = QSM_simple_clean.cylinder.start;
        plot3(start_clean(:,1), start_clean(:,2), start_clean(:,3), '.k', 'MarkerSize', 4);
    end
    xlabel('x [m]'); ylabel('y [m]'); zlabel('z [m]');
    title('Simplified model (no islands)');
    axis equal
    grid off
    view(3)
end

%% 17) VOLUME TABLE with variability, in m^3------------------------
%  
%       All inputs              = all models of the first run
%       Optimal                 = models of the winning parameter combination
%       Optimal (single)        = the one model that gets simplified
%       Estimated               = optimal models + second run (better std)
%       Simplified              = the simplified model (single -> no std)
%       Simplified (no islands) = Simplified, with disconnected branch
%                                  islands removed (section 15c) - only
%                                  present if islands were found
%       TreeQSM stores volumes in LITERS -> divide by 1000
%  ------------------------------------------------------------

% Takes an array of QSMs, returns an n-by-3 matrix [total, stem, branch] in m^3
vols = @(QA) cell2mat(arrayfun(@(Q) ...
    [Q.treedata.TotalVolume, Q.treedata.TrunkVolume, Q.treedata.BranchVolume] ./ 1000, ...
    QA(:), 'UniformOutput', false));

% --- volume after filtering out cylinders thinner than cut_cm ----
% (moved here from former section 17b, so the result can go into VolumeTable)
cut_cm = 10;                      % cut-off diameter [cm]

r_cyl = QSM_opt.cylinder.radius;      % [m]
L_cyl = QSM_opt.cylinder.length;      % [m]
V_cyl = pi .* r_cyl.^2 .* L_cyl;      % [m^3]
d_cm  = 2 * r_cyl * 100;              % diameter [cm]

keep = d_cm >= cut_cm;

fprintf('\n--- Cylinders, cut-off %.0f cm ---\n', cut_cm);
fprintf('Cylinders total  : %d\n', numel(r_cyl));
fprintf('Cylinders kept   : %d (%.1f %%)\n', sum(keep), sum(keep)/numel(r_cyl)*100);
fprintf('Volume total     : %.3f m3\n', sum(V_cyl));
fprintf('Volume kept      : %.3f m3\n', sum(V_cyl(keep)));
fprintf('Volume removed   : %.3f m3 (%.1f %%)\n', ...
    sum(V_cyl(~keep)), sum(V_cyl(~keep))/sum(V_cyl)*100);

% --- split cylinders into stem (order 0) and branches (order >= 1) ---
order = QSM_opt.cylinder.BranchOrder;   % 0 = stem, >=1 = branch

is_stem   = (order == 0);
is_branch = (order >= 1);

% combine with the diameter filter (keep) computed above
keep_stem   = keep & is_stem;
keep_branch = keep & is_branch;

Vstem_filt   = sum(V_cyl(keep_stem));
Vbranch_filt = sum(V_cyl(keep_branch));
Vtotal_filt  = sum(V_cyl(keep));     % should equal Vstem_filt + Vbranch_filt

fprintf('Stem volume kept    : %.3f m3\n', Vstem_filt);
fprintf('Branch volume kept  : %.3f m3\n', Vbranch_filt);

% --- unfiltered stem/branch volumes for comparison ---
Vstem_total   = sum(V_cyl(is_stem));
Vbranch_total = sum(V_cyl(is_branch));

Vstem_removed_pct   = (Vstem_total   - Vstem_filt)   / Vstem_total   * 100;
Vbranch_removed_pct = (Vbranch_total - Vbranch_filt) / Vbranch_total * 100;

fprintf('Stem volume removed  : %.3f m3 (%.1f %%)\n', ...
    Vstem_total - Vstem_filt, Vstem_removed_pct);
fprintf('Branch volume removed: %.3f m3 (%.1f %%)\n', ...
    Vbranch_total - Vbranch_filt, Vbranch_removed_pct);

% --- same cut-off, but applied to LENGTH instead of volume ------------
% Uses the SAME keep/is_stem/is_branch masks and the SAME L_cyl array
% already used for Vstem_filt/Vbranch_filt above - just sum(L_cyl(...))
% instead of sum(V_cyl(...)). This lets "TreeQSM mine (*, Filtered<10cm)"
% report its OWN filtered trunk/branch length (exported in section 18
% below) instead of reusing the unfiltered model's length.
Lstem_filt   = sum(L_cyl(keep_stem));
Lbranch_filt = sum(L_cyl(keep_branch));

Lstem_total   = sum(L_cyl(is_stem));
Lbranch_total = sum(L_cyl(is_branch));

Lstem_removed_pct   = (Lstem_total   - Lstem_filt)   / Lstem_total   * 100;
Lbranch_removed_pct = (Lbranch_total - Lbranch_filt) / Lbranch_total * 100;

fprintf('Stem length kept     : %.3f m\n', Lstem_filt);
fprintf('Branch length kept   : %.3f m\n', Lbranch_filt);
fprintf('Stem length removed  : %.3f m (%.1f %%)\n', ...
    Lstem_total - Lstem_filt, Lstem_removed_pct);
fprintf('Branch length removed: %.3f m (%.1f %%)\n', ...
    Lbranch_total - Lbranch_filt, Lbranch_removed_pct);

V_filtered = [Vtotal_filt, Vstem_filt, Vbranch_filt];
fprintf('Check: stem+branch total = %.3f m3 (should equal %.3f m3)\n', ...
    Vstem_total + Vbranch_total, sum(V_cyl));
%
% ncyl(QA): companion to vols() above - returns an n-by-1 vector with the
% number of cylinders in each model of QA (one count per model, same row
% order as vols(QA)). Used to build the new "n_cylinders" column of
% VolumeTable below (Task C): every model in a group gets its cylinder
% count averaged into that group's row, the same way vols() values get
% averaged into Mean_m3 - so cylinder count is tracked "for all models",
% not just a single representative one.
ncyl = @(QA) arrayfun(@(Q) numel(Q.cylinder.radius), QA(:));

% --- collect the groups --------------------------------------
V_all  = vols(res.QSMs);              % whole parameter grid
V_opt  = vols(res.QSMs(idx));         % winning combination
V_rep  = vols(res.QSMs(idx_rep));     % the single model that gets simplified
V_simp = vols(QSM_simple(end));       % simplified model

% Ncyl_* vectors mirror V_* above, row-for-row (same models, same order) -
% this is what lets the table-build loop just do mean(Ncyl_*) the same
% way it already does mean(V_*).
Ncyl_all = ncyl(res.QSMs);
Ncyl_opt = ncyl(res.QSMs(idx));
Ncyl_rep = ncyl(res.QSMs(idx_rep));

% groups now has a THIRD column: the per-model cylinder-count vector
% matching that row's V_matrix (Task C). Column order stays {name, V, Ncyl}.
groups = {'All inputs',       V_all,  Ncyl_all;
          'Optimal',          V_opt,  Ncyl_opt;
          'Optimal (single)', V_rep,  Ncyl_rep};

% Estimated = optimal group + second run, for a better std estimate
if exist('res_new', 'var')
    V_est = [V_opt; vols(res_new.QSMs)];
    Ncyl_est = [Ncyl_opt; ncyl(res_new.QSMs)];   % same concatenation as V_est, so rows still line up
    groups(end+1,:) = {'Estimated', V_est, Ncyl_est};
    fprintf('Estimated group: %d models.\n', size(V_est,1));
else
    warning('res_new not found - run steps 13 and 14 to get the Estimated row.');
end

% Keep 'Simplified' (uncleaned QSM_simple) as-is for side-by-side
% comparison against the cleaned version added right below.
Ncyl_simp = ncyl(QSM_simple(end));
groups(end+1,:) = {'Simplified', V_simp, Ncyl_simp};

% QSM_simple_clean only exists if section 15c actually removed islands
% (it is NOT run/created when island_groups was empty - see its comment
% above). This adds a SEPARATE new row ('Simplified (no islands)')
% without touching or replacing the existing 'Simplified' row above, so
% both are kept side by side in VolumeTable/volumes_*.csv for comparison.
if exist('QSM_simple_clean', 'var')
    V_simp_clean = vols(QSM_simple_clean(end));   % same vols() helper used for V_simp above
    Ncyl_simp_clean = ncyl(QSM_simple_clean(end));
    groups(end+1,:) = {'Simplified (no islands)', V_simp_clean, Ncyl_simp_clean};
else
    fprintf('QSM_simple_clean not found - skipping "Simplified (no islands)" group (no islands detected, or section not run).\n');
end

% Only add the "Filtered <10cm" group when export_filtered_10cm is true
% (set in section 1c). This is the switch itself: if it's false, this row
% is simply never appended to `groups`, so it never reaches VolumeTable.
if export_filtered_10cm
    % This group has no separate QSM struct of its own - it's just QSM_opt's
    % cylinders restricted by the "keep" mask (built earlier in this same
    % section, diameter >= cut_cm). So its cylinder count is simply how
    % many entries in "keep" are true, i.e. sum(keep) - already the exact
    % same count used for Cylinders kept/Vtotal_filt above.
    Ncyl_filtered = sum(keep);
    groups(end+1,:) = {'Filtered <10cm', V_filtered, Ncyl_filtered};
end

% --- build the table -----------------------------------------
attr = ["Total"; "Stem"; "Branches"];

Group = strings(0,1); Attribute = strings(0,1);
N = []; Mean_m3 = []; Std_m3 = []; CV_pct = [];
N_cylinders = [];   % Task C: new column, average cylinder count for this group

for g = 1:size(groups,1)
    V = groups{g,2};
    n = size(V,1);                 % number of models in this group
    m = mean(V, 1);                % mean of each column
    if n > 1
        s = std(V, 0, 1);          % sample standard deviation
    else
        s = nan(1,3);              % single model -> std undefined
    end
    cv = s ./ m .* 100;            % coefficient of variation [%]

    % Average cylinder count across every model in this group, same idea
    % as mean(V,1) above but for Ncyl - round() because a cylinder count
    % is always a whole number, and averaging several models' counts can
    % otherwise land on a fraction (e.g. 5 models with 100/101 cylinders).
    ncyl_mean = round(mean(groups{g,3}));

    for a = 1:3
        Group(end+1,1)       = groups{g,1};
        Attribute(end+1,1)   = attr(a);
        N(end+1,1)           = n;
        Mean_m3(end+1,1)     = m(a);
        Std_m3(end+1,1)      = s(a);
        CV_pct(end+1,1)      = cv(a);
        % Same value repeated for all 3 attribute rows of this group - a
        % cylinder count belongs to the whole model, not to one specific
        % attribute (Total/Stem/Branches), same pattern already used for N.
        N_cylinders(end+1,1) = ncyl_mean;
    end
end

Tree = repmat(string(tree_id), height(Group), 1);
Run  = repmat(string(run_tag), height(Group), 1);
% N_cylinders added as the LAST column (Task C) - existing columns/order
% (Tree..CV_pct) are unchanged so any code still expecting the old shape
% keeps working, it just also gets this one extra column now.
VolumeTable = table(Tree, Run, Group, Attribute, N, Mean_m3, Std_m3, CV_pct, N_cylinders);
VolumeTable.Properties.Description = tree_id;
disp(VolumeTable);

% --- reference values of the optimal model -------------------
fprintf('DBHqsm     = %.1f cm\n', QSM_opt.treedata.DBHqsm * 100);
fprintf('DBHcyl     = %.1f cm\n', QSM_opt.treedata.DBHcyl * 100);
fprintf('TreeHeight = %.2f m\n',  QSM_opt.treedata.TreeHeight);

% --- save ----------------------------------------------------
save(vol_file, 'VolumeTable');
writetable(VolumeTable, [vol_file '.csv']);
fprintf('Volume table saved as %s.mat and %s.csv\n', vol_file, vol_file);
%% 18) DBH AND HEIGHT FOR COMPARISON----------------
%  
%  ------------------------------------------------------------
dbh_file = ['dbh_' tree_id '_' run_tag '.txt'];
fid = fopen(dbh_file, 'w');
fprintf(fid, '%.6f', QSM_opt.treedata.DBHqsm);   % stem diameter at 1.3 m [m]
fclose(fid);
fprintf('DBH exported to %s\n', dbh_file);
%
% ---- export tree height for the shared results table ----
h_file = ['height_' tree_id '_' run_tag '.txt'];
fid = fopen(h_file, 'w');
fprintf(fid, '%.6f', QSM_opt.treedata.TreeHeight);   % tree height [m]
fprintf('Height exported to %s\n', h_file);
fclose(fid);
%
% ---- export TAPER (stem narrowing) for the shared results table ----
% taper_cm_per_m = (diameter_at_1.3m - diameter_at_10.0m) * 100 / (10.0 - 1.3)
%
% 1.3 m and 10.0 m are NOT arbitrary here - they are exactly TAPER_H_LOWER/
% TAPER_H_UPPER, the same two reference heights every OTHER method in this
% pipeline uses for its own taper_cm_per_m (see adtree_reconstruct_compare.py's
% raw_taper/cal_taper, and TAPER_H_LOWER/TAPER_H_UPPER in
% import_matlab_results.py) - using the same two heights here is what makes
% this number directly comparable to theirs in volume_results.csv, instead
% of being a taper measured over a different (and therefore not comparable)
% span of the trunk.
%
% Diameter at each height is found the SAME way QSM_opt.treedata.DBHqsm
% itself is computed (see dbh_cylinder() in main_steps/tree_data.m): walk
% the TRUNK cylinders in order from the base upward, add up their lengths,
% and take the diameter of the first cylinder whose running total reaches
% the target height. That original code walks the raw trunk point cloud
% (no longer available here, this far into the script) - but section 17b
% above already has the exact same trunk cylinder radii/lengths in scope
% (r_cyl, L_cyl, is_stem), in the same base-to-tip order, so we can repeat
% the identical "which cylinder is at height h" logic using THOSE arrays
% instead, just evaluated at two heights (1.3 m and 10.0 m) rather than one.
trunk_len_cyl = L_cyl(is_stem);     % trunk-only cylinder lengths, base -> tip order
trunk_rad_cyl = r_cyl(is_stem);     % trunk-only cylinder radii, same order
cum_h = cumsum(trunk_len_cyl);      % running height from the base, one value per trunk cylinder

idx_13 = find(cum_h >= 1.3,  1, 'first');   % index of the first trunk cylinder reaching 1.3 m
idx_10 = find(cum_h >= 10.0, 1, 'first');   % index of the first trunk cylinder reaching 10.0 m

taper_file = ['taper_' tree_id '_' run_tag '.txt'];
fid = fopen(taper_file, 'w');
if ~isempty(idx_13) && ~isempty(idx_10)
    d_13 = 2 * trunk_rad_cyl(idx_13);   % stem diameter [m] at 1.3 m
    d_10 = 2 * trunk_rad_cyl(idx_10);   % stem diameter [m] at 10.0 m
    taper_cm_per_m = (d_13 - d_10) * 100.0 / (10.0 - 1.3);
    fprintf(fid, '%.6f', taper_cm_per_m);
    fprintf('Taper exported to %s (%.2f cm/m)\n', taper_file, taper_cm_per_m);
    % Sanity check only (not written to file): d_13 should equal
    % QSM_opt.treedata.DBHqsm, since it's the exact same walking algorithm
    % applied to the exact same trunk cylinders - if these two numbers
    % don't (nearly) match when you run this, the array order assumption
    % above is wrong for this tree and the taper value should not be trusted.
    fprintf('  (check: diameter at 1.3 m = %.4f m, vs. treedata.DBHqsm = %.4f m)\n', ...
        d_13, QSM_opt.treedata.DBHqsm);
else
    % The trunk model doesn't reach 10.0 m (short tree, or the trunk
    % cylinders stop earlier than that) - write an EMPTY file rather than
    % extrapolating or guessing a number. import_matlab_results.py's
    % read_single_number() already treats an empty/unparseable file as
    % "value not available" (None), exactly like every other optional
    % field in this pipeline (DBH/height/trunk-length/branch-length).
    fprintf(fid, '');
    fprintf('Taper NOT exported - trunk model does not reach 10.0 m (empty %s written)\n', taper_file);
end
fclose(fid);
%
% ---- export trunk/branch length for the shared results table ----
% QSM_opt.treedata.TrunkLength/BranchLength are set in tree_data.m as
% sum(Len(Trunk)) / sum(Len(~Trunk)) - i.e. the exact same cylinder length
% array used to compute TrunkVolume/BranchVolume above, just summed instead
% of pi*r^2*length-weighted-summed. Already in metres (unlike the *Volume
% fields, which are in litres), so no unit conversion is needed here.
trunklen_file = ['trunklen_' tree_id '_' run_tag '.txt'];
fid = fopen(trunklen_file, 'w');
fprintf(fid, '%.6f', QSM_opt.treedata.TrunkLength);   % trunk/stem length [m]
fprintf('Trunk length exported to %s\n', trunklen_file);
fclose(fid);
%
branchlen_file = ['branchlen_' tree_id '_' run_tag '.txt'];
fid = fopen(branchlen_file, 'w');
fprintf(fid, '%.6f', QSM_opt.treedata.BranchLength);   % branch length [m]
fprintf('Branch length exported to %s\n', branchlen_file);
fclose(fid);
%
% ---- export FILTERED (>= cut_cm diameter) trunk/branch length ----
% ADDED, does not replace trunklen_*.txt/branchlen_*.txt above: those two
% still hold the UNFILTERED length (used by Optimal/Estimated/... groups).
% These two new files hold Lstem_filt/Lbranch_filt from section 17b, for
% the "Filtered <10cm" group specifically - different filename so neither
% overwrites the other, and import_matlab_results.py picks the right one
% based on which group it's importing.
% Same switch as the "Filtered <10cm" VolumeTable row above (section 1c):
% when export_filtered_10cm is false, skip writing these two text files
% entirely, instead of writing files nobody asked for on this run.
if export_filtered_10cm
    trunklen_filtered_file = ['trunklen_filtered_' tree_id '_' run_tag '.txt'];
    fid = fopen(trunklen_filtered_file, 'w');
    fprintf(fid, '%.6f', Lstem_filt);   % filtered trunk/stem length [m]
    fprintf('Filtered trunk length exported to %s\n', trunklen_filtered_file);
    fclose(fid);
    %
    branchlen_filtered_file = ['branchlen_filtered_' tree_id '_' run_tag '.txt'];
    fid = fopen(branchlen_filtered_file, 'w');
    fprintf(fid, '%.6f', Lbranch_filt);   % filtered branch length [m]
    fprintf('Filtered branch length exported to %s\n', branchlen_filtered_file);
    fclose(fid);
else
    fprintf('export_filtered_10cm = false - skipping trunklen_filtered_*.txt and branchlen_filtered_*.txt\n');
end
%% 19) EXPORT - the actual reconstruction parameters used for this run ----
%  ------------------------------------------------------------
% Sidecar CSV (header + ONE data row) recording the ACTUAL values used for
% THIS run - inputs.PatchDiam1/PatchDiam2Min/PatchDiam2Max/MinCylRad hold
% the right value regardless of whether manual_patchdiam was true or false:
% the auto branch (section 9's "else") never touches inputs.PatchDiam1/
% PatchDiam2Min/PatchDiam2Max/MinCylRad at all, so they're always readable
% straight off `inputs` here either way (see section 9's "define_input
% (automatic)" vs. "MANUAL PatchDiam" blocks above) - no special-casing
% needed for the auto case.
%
% simp_MaxOrder/simp_SmallRadii/simp_ReplaceIterations are always set
% (section 16), regardless of manual/auto mode, so they're read straight
% from those variables rather than from `inputs`.
%
% tree/run use the SAME tree_id/run_tag values already used to build
% vol_file/export_prefix (section 2) - the whole point of this file is to
% let import_matlab_results.py join it back to the right volumes_*.csv row
% by that same (tree, run) pair, read ONCE per file (these values don't
% vary by Group the way total/trunk/branch volume does - they were fixed
% for the whole MATLAB run).
mode_str = 'auto';
if manual_patchdiam
    mode_str = 'manual';
end
params_file = ['params_' tree_id '_' run_tag '.csv'];
fid = fopen(params_file, 'w');
fprintf(fid, 'tree,run,mode,pd1_m,pd2min_m,pd2max_m,mincylrad_m,simp_maxorder,simp_smallradii,simp_replaceiterations\n');
fprintf(fid, '%s,%s,%s,%.6f,%.6f,%.6f,%.6f,%d,%.6f,%d\n', ...
    tree_id, run_tag, mode_str, ...
    inputs.PatchDiam1, inputs.PatchDiam2Min, inputs.PatchDiam2Max, inputs.MinCylRad, ...
    simp_MaxOrder, simp_SmallRadii, simp_ReplaceIterations);
fclose(fid);
fprintf('Parameters exported to %s\n', params_file);

%% 20) EXPORT geometry for ANSYS----------------
%  ------------------------------------------------------------
% Gated behind export_to_ansys (section 1c, default false) - same
% "explicit opt-in" pattern as clean_start/export_filtered_10cm - so this
% does NOT run on every single cell-by-cell iteration, only when you
% actually want geometry exported for this run.
if export_to_ansys

% EXPORT_FROM_SAVED_RUN_TAG (section 1c): non-blank means "load a
% previously-saved simplified_<tree_id>_<tag>.mat from disk (section 16c's
% unconditional save) and export THAT" - skips the ansys_source switch
% below entirely, ignoring whatever's currently live in memory. Blank
% (default) keeps today's behaviour: export whatever's live in memory,
% via the existing ansys_source switch.
if ~isempty(EXPORT_FROM_SAVED_RUN_TAG)
    saved_file = ['simplified_' tree_id '_' EXPORT_FROM_SAVED_RUN_TAG '.mat'];
    if ~isfile(saved_file)
        error('EXPORT_FROM_SAVED_RUN_TAG = ''%s'' but %s does not exist.', ...
              EXPORT_FROM_SAVED_RUN_TAG, saved_file);
    end
    loaded = load(saved_file);
    qsm_selected = loaded.QSM_final;
    n_opt = 1;
    ansys_tag_default = EXPORT_FROM_SAVED_RUN_TAG;
    fprintf('Exporting from SAVED variant: %s\n', saved_file);
else
    % --- choose the SOURCE model ------------------------------------
    % 'simplified_clean' = QSM_simple_clean (simplified model, island cylinders
    %                       removed - section 15c). This is now the PRIMARY
    %                       model the user exports to ANSYS: it is both the
    %                       fewer-cylinder simplified geometry AND has any
    %                       disconnected branch-island noise stripped out.
    % 'simplified'        = QSM_simple (after simplify_qsm, section 16, BEFORE
    %                        island cleaning) - kept for occasional comparison
    %                        exports only.
    % 'optimal'           = QSM_opt (before simplification, section 16, output
    %                        of select_optimum) - kept for occasional
    %                        comparison exports only.
    ansys_source = 'simplified_clean';   % <-- default: simplified + islands removed

    switch ansys_source
        case 'simplified_clean'
            % QSM_simple_clean only exists if section 15c actually ran AND
            % found islands to remove (see its comment above). If it wasn't
            % created, error out with a clear pointer to 'simplified' instead
            % of silently exporting the wrong thing - in the no-islands case
            % QSM_simple and QSM_simple_clean would be identical anyway, so
            % nothing is lost by switching ansys_source manually.
            if ~exist('QSM_simple_clean', 'var')
                error(['ansys_source = ''simplified_clean'', but QSM_simple_clean does not exist ' ...
                       '(no islands were found in section 15b/15c, so there was nothing to clean). ' ...
                       'Set ansys_source = ''simplified'' instead - QSM_simple and QSM_simple_clean ' ...
                       'would be identical anyway when no islands exist.']);
            end
            qsm_selected = QSM_simple_clean(1);   % single cleaned model (not an array)
        case 'simplified'
            ansys_export_idx = 1;                 % index into QSM_simple
            qsm_selected = QSM_simple(ansys_export_idx);
        case 'optimal'
            qsm_selected = QSM_opt;               % QSM_opt is a single model (not an array)
        otherwise
            error('ansys_source must be ''simplified_clean'', ''simplified'' or ''optimal''.');
    end

    n_opt = length(qsm_selected);

    fprintf('Exporting to ANSYS from source: %s (%d model(s))\n', ansys_source, n_opt);
    ansys_tag_default = run_tag;
end

% ansys_tag: ansys_export_name (section 1c) overrides ansys_tag_default
% (EXPORT_FROM_SAVED_RUN_TAG when exporting a saved variant, run_tag
% otherwise) when non-blank, so a one-off export can be clearly labeled
% with something more memorable; blank (the default) just falls back to
% ansys_tag_default, same as every other export in this file.
if isempty(ansys_export_name)
    ansys_tag = ansys_tag_default;
else
    ansys_tag = ansys_export_name;
end
% Built the SAME way export_prefix is built in section 2, just from
% ansys_tag instead of run_tag directly - so a custom ansys_export_name
% actually changes the exported filenames, not just a label.
ansys_export_prefix = ['geom_' tree_id '_' ansys_tag '_'];

geom_orig = myfun.result_ansys(qsm_selected, n_opt);

for i = 1:n_opt
    geom_table = geom_orig{i};                                   % table of one model
    file_name  = sprintf('%s%d.txt', ansys_export_prefix, i);   % e.g. geom_IND07_v3_1.txt
    writematrix(geom_table, file_name, 'Delimiter', '\t');
    fprintf('Exported: %s\n', file_name);
end

else
    fprintf('export_to_ansys = false - skipping ANSYS geometry export.\n');
end

% REMOVED: a duplicate, leftover second "20) EXPORT geometry for ANSYS"
% block used to sit here. It always re-exported QSM_simple (ignoring
% whatever ansys_source above had chosen, e.g. 'optimal') using
% n_opt = length('QSM_simple') - a bug, since that computes the length of
% the literal 10-character STRING 'QSM_simple', not length(QSM_simple)
% (the number of models). Because it wrote to the exact same file_name
% pattern as the block above, it silently ran a second time and could
% overwrite the export the user actually intended. Deleted as dead code -
% the switch above already covers every case correctly on its own.

% ---------------------------------------------------------------
% LOCAL FUNCTIONS
% A script file (as opposed to a function file) is allowed to define its
% own "local functions" (small helper functions usable only inside this
% same file), but ONLY if they are placed at the very end of the file,
% after every other script statement - MATLAB would otherwise not know
% where the script code ends and the function definitions begin.
% ---------------------------------------------------------------

function tag = compute_run_tag(manual_patchdiam, man_PD1, man_PD2Min, ...
                                 man_PD2Max, simp_MaxOrder, simp_SmallRadii, ...
                                 simp_ReplaceIterations)
    % The ONE place run_tag's formula lives - called from section 2 (to
    % actually SET run_tag) and from section 16's sanity check (to
    % recompute what run_tag SHOULD be from current settings and compare)
    % - so the two can never drift out of sync the way two independently
    % hand-copied formulas eventually would.
    if manual_patchdiam
        mode_tag = sprintf('man_pd%02d-%02d-%02d', round(man_PD1*100), ...
            round(man_PD2Min*100), round(man_PD2Max*100));
    else
        mode_tag = 'aut';
    end
    simp_tag = sprintf('mo%d_sr%03d_ri%d', simp_MaxOrder, ...
        round(simp_SmallRadii*1000), simp_ReplaceIterations);
    tag = [mode_tag '_' simp_tag];
end

function island_groups = find_disconnected_islands(parent_arr)
    % Finds every group of cylinders that got disconnected from the
    % main tree structure (cylinder.parent == 0 for any cylinder OTHER
    % than cylinder index 1, which is always the true trunk-base root).
    % Returns a cell array, one entry per island, each containing the
    % list of cylinder indices belonging to that island (the orphan
    % root itself plus every descendant found by walking DOWN the
    % parent-child tree from it).
    n_cyl = numel(parent_arr);
    children = cell(n_cyl, 1);
    for k = 1:n_cyl
        p = parent_arr(k);
        if p > 0
            children{p} = [children{p}, k];
        end
    end

    orphan_roots = find(parent_arr == 0);
    orphan_roots = orphan_roots(orphan_roots ~= 1);   % exclude the true root

    island_groups = cell(numel(orphan_roots), 1);
    for oi = 1:numel(orphan_roots)
        stack = orphan_roots(oi);
        island_indices = [];
        while ~isempty(stack)
            cur = stack(end);
            stack(end) = [];
            island_indices(end+1) = cur;
            stack = [stack, children{cur}];
        end
        island_groups{oi} = island_indices;
    end
end