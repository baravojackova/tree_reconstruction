%% ============================================================
%  TreeQSM - reconstruction of tree structure from a point cloud
%  All file names are defined in ONE place (USER SETTINGS below).
%  Run the script cell by cell (Ctrl+Enter in each %% section).
%  ============================================================
clear
clc

%% ------------------------------------------------------------
%  1) USER SETTINGS - this is the ONLY block you need to edit
%  ------------------------------------------------------------
run_tag = 'v2manual';    % change for every new settings variant

% Folder with the TreeQSM source code (contains treeqsm.m, +myfun, ...)
src_dir = 'C:\Users\Spravce\Documents\BARA\01_Skeny_Babice\tree_reconstruction\scripts\matlab';

% Folder with the point clouds and where all results will be written
data_dir = 'C:\Users\Spravce\Documents\BARA\01_Skeny_Babice\tree_reconstruction\scripts\matlab';

% --- tree identification -------------------------------------
tree_id   = 'IND01_054';           % short name used for ALL output files
cloud_txt = 'IND01_054.txt';   % input point cloud (text file, 3 columns X Y Z)

% --- number of models ----------------------------------------
n_models_first = 5;    % models per parameter combination, first (coarse) run
n_models_opt   = 25;   % models with the optimal inputs, second run

% --- parallel computing --------------------------------------
use_parallel = true;   % true  = make_models_parallel
                       % false = make_models (single core, fallback)
n_workers = 0;         % 0 = derive automatically from the number of tasks

% --- parameter ranges for define_input -----------------------
% define_input(P, nPD1, nPD2Min, nPD2Max) = how many values are tested
% for each of the three patch-diameter parameters
% !!! if manual input nPD = 1
nPD1    = 1;
nPD2Min = 1;
nPD2Max = 1;

% --- MANUAL PatchDiam (see section 9) ------------------------
manual_patchdiam = true;   % false = keep everything from define_input

% PatchDiam1 (rought first cover) has to be t ≥ PatchDiam2Max (gentle cover)

man_PD1    = 0.08;   % PatchDiam1     - AdQSM paper, Indonesian site 0,08
man_PD2Min = 0.02;   % PatchDiam2Min
man_PD2Max = 0.07;   % PatchDiam2Max
man_MinCylRad = 0.0025; % MinCylRad

% --- which model is simplified and exported ------------------
use_optimal = true;    % true  = use the optimal model from select_optimum
                       % false = use res.QSMs(model_index)
model_index = 1;       % used only when use_optimal = false

plot_optimal = true;   % true = plot the optimal QSM before simplification

% --- simplification settings ---------------------------------
simp_MaxOrder          = 10;
simp_SmallRadii        = 0.05;
simp_ReplaceIterations = 1;
simp_Plot              = 1;
simp_Disp              = 1;

%% ------------------------------------------------------------
%  2) DERIVED NAMES - built automatically from tree_id + run_tag
%     Normally you do NOT edit this block
%  ------------------------------------------------------------
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

%% ------------------------------------------------------------
%  3) IMPORT the point cloud from the text file
%  ------------------------------------------------------------
if ~isfile(cloud_txt)
    error('Input file not found: %s', fullfile(data_dir, cloud_txt));
end
P = load(cloud_txt);
fprintf('Loaded %d points.\n', size(P,1));

%% ------------------------------------------------------------
%  4) SHIFT the coordinate system so that the cloud starts at 0
%  ------------------------------------------------------------
P(:,1) = P(:,1) - min(P(:,1));   % X
P(:,2) = P(:,2) - min(P(:,2));   % Y
P(:,3) = P(:,3) - min(P(:,3));   % Z

fprintf('Extent: X %.2f m, Y %.2f m, Z %.2f m\n', ...
    max(P(:,1)), max(P(:,2)), max(P(:,3)));

%% ------------------------------------------------------------
%  5) PLOT the point cloud for visual check
%  ------------------------------------------------------------
figure('Name', ['check - ' tree_id]);
plot3(P(:,1), P(:,2), P(:,3), '.k', 'MarkerSize', 1);
axis equal;                      % same scale on all axes
grid on;
xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]');
title(['Cleaned and positioned tree: ' tree_id], 'Interpreter', 'none');

%% ------------------------------------------------------------
%  6) SAVE the point cloud as a .mat file
%  ------------------------------------------------------------
save(mat_file, 'P');
fprintf('Point cloud saved as %s\n', mat_file);

%% ------------------------------------------------------------
%  7) RE-ENTRY POINT - start here if the .mat file already exists
%  ------------------------------------------------------------
load(mat_file);                  % loads variable P

%% ------------------------------------------------------------
%  8) DEFINE the input parameters automatically
%     define_input first runs create_input, so ALL other fields
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

%% ------------------------------------------------------------
%  9) MANUAL PatchDiam - BallRad derived with the SAME formulas
%     that define_input uses:
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

%% ------------------------------------------------------------
%  9b) START the parallel pool (run this BEFORE step 10)
%      Doing it here separates pool problems from QSM problems
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

%% ------------------------------------------------------------
%  10) FIRST RUN - models over the parameter grid
%  ------------------------------------------------------------
tic     % start time monitoring
if use_parallel
    QSMs = make_models_parallel(mat_name, res_name, n_models_first, inputs);
else
    QSMs = make_models(mat_name, res_name, n_models_first, inputs);
end
fprintf('First QSM run finished in %.1f min.\n', toc/60);   % toc - stop time moniroting

%% ------------------------------------------------------------
%  11) LOAD the results of the first run
%  ------------------------------------------------------------
res = load(res_file);
fprintf('Loaded %d models from %s\n', numel(res.QSMs), res_file);

%% ------------------------------------------------------------
%  12) OPTIMISATION - select the best models
%      (based on point-to-cylinder distances)
%  ------------------------------------------------------------
[TreeData_O, OptModels, OptInputs, OptQSM] = select_optimum(res.QSMs);

%% ------------------------------------------------------------
%  13) SECOND RUN - more models with the optimal parameters
%  ------------------------------------------------------------
tic
if use_parallel
    QSMs_new = make_models_parallel(mat_name, res_new_name, n_models_opt, OptInputs);
else
    QSMs_new = make_models(mat_name, res_new_name, n_models_opt, OptInputs);
end
fprintf('Second QSM run finished in %.1f min.\n', toc/60);

%% ------------------------------------------------------------
%  14) LOAD the results of the second run
%  ------------------------------------------------------------
res_new = load(res_new_file);
fprintf('Loaded %d models from %s\n', numel(res_new.QSMs), res_new_file);

%% ------------------------------------------------------------
%  15) PRECISION - combine both runs to get better std estimates
%  ------------------------------------------------------------
[TreeData_E, OptQSMs_E, OptQSM_E] = estimate_precision( ...
    res.QSMs, res_new.QSMs, TreeData_O, OptModels);

%% ------------------------------------------------------------
%  16) SELECT the source model, then SIMPLIFY it
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

QSM_simple = simplify_qsm(QSM_opt, simp_MaxOrder, ...
    simp_SmallRadii, simp_ReplaceIterations, simp_Plot, simp_Disp);

%% ------------------------------------------------------------
%  17) VOLUME TABLE with variability, in m^3
%       All inputs       = all models of the first run
%       Optimal          = models of the winning parameter combination
%       Optimal (single) = the one model that gets simplified
%       Estimated        = optimal models + second run (better std)
%       Simplified       = the simplified model (single -> no std)
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
% --- collect the groups --------------------------------------
V_all  = vols(res.QSMs);              % whole parameter grid
V_opt  = vols(res.QSMs(idx));         % winning combination
V_rep  = vols(res.QSMs(idx_rep));     % the single model that gets simplified
V_simp = vols(QSM_simple(end));       % simplified model

groups = {'All inputs',       V_all;
          'Optimal',          V_opt;
          'Optimal (single)', V_rep};

% Estimated = optimal group + second run, for a better std estimate
if exist('res_new', 'var')
    V_est = [V_opt; vols(res_new.QSMs)];
    groups(end+1,:) = {'Estimated', V_est};
    fprintf('Estimated group: %d models.\n', size(V_est,1));
else
    warning('res_new not found - run steps 13 and 14 to get the Estimated row.');
end

groups(end+1,:) = {'Simplified', V_simp};
groups(end+1,:) = {'Filtered <10cm', V_filtered};

% --- build the table -----------------------------------------
attr = ["Total"; "Stem"; "Branches"];

Group = strings(0,1); Attribute = strings(0,1);
N = []; Mean_m3 = []; Std_m3 = []; CV_pct = [];

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

    for a = 1:3
        Group(end+1,1)     = groups{g,1};
        Attribute(end+1,1) = attr(a);
        N(end+1,1)         = n;
        Mean_m3(end+1,1)   = m(a);
        Std_m3(end+1,1)    = s(a);
        CV_pct(end+1,1)    = cv(a);
    end
end

Tree = repmat(string(tree_id), height(Group), 1);
Run  = repmat(string(run_tag), height(Group), 1);
VolumeTable = table(Tree, Run, Group, Attribute, N, Mean_m3, Std_m3, CV_pct);
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
%% ------------------------------------------------------------
%  18) DBH AND HEIGHT FOR COMPARISON
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
%% ------------------------------------------------------------
%  19) EXPORT - table with the input parameters of each model
%  ------------------------------------------------------------
n = length(QSM_simple);
info = table('Size', [n 7], ...
    'VariableTypes', {'string','double','double','double','double','double','double'}, ...
    'VariableNames', {'Name','tree','model','PD1','PD2Min','PD2Max','Time'});

for i = 1:n
    info.Name(i)   = QSM_simple(i).rundata.inputs.name;
    info.tree(i)   = QSM_simple(i).rundata.inputs.tree;
    info.model(i)  = QSM_simple(i).rundata.inputs.model;
    info.PD1(i)    = QSM_simple(i).rundata.inputs.PatchDiam1;
    info.PD2Min(i) = QSM_simple(i).rundata.inputs.PatchDiam2Min;
    info.PD2Max(i) = QSM_simple(i).rundata.inputs.PatchDiam2Max;
    info.Time(i)   = QSM_simple(i).rundata.time(end,1)./60;   % minutes
end

time_sum      = sum(info.Time);       % total time in minutes
time_sum(1,2) = time_sum./60;         % total time in hours

%% ------------------------------------------------------------
%  20) EXPORT geometry for ANSYS
%  ------------------------------------------------------------
% --- choose the SOURCE model ------------------------------------
% 'simplified' = QSM_simple (after simplify_qsm, section 16)
% 'optimal'    = QSM_opt    (before simplification, section 16, output of select_optimum)
ansys_source = 'optimal';   % <-- change to 'simplified' when you want the simplified export

switch ansys_source
    case 'simplified'
        ansys_export_idx = 1;                 % index into QSM_simple
        qsm_selected = QSM_simple(ansys_export_idx);
    case 'optimal'
        qsm_selected = QSM_opt;               % QSM_opt is a single model (not an array)
    otherwise
        error('ansys_source must be ''simplified'' or ''optimal''.');
end

n_opt = length(qsm_selected);

fprintf('Exporting to ANSYS from source: %s (%d model(s))\n', ansys_source, n_opt);

geom_orig = myfun.result_ansys(qsm_selected, n_opt);

for i = 1:n_opt
    geom_table = geom_orig{i};                            % table of one model
    file_name  = sprintf('%s%d.txt', export_prefix, i);   % e.g. geom_IND07_v3_1.txt
    writematrix(geom_table, file_name, 'Delimiter', '\t');
    fprintf('Exported: %s\n', file_name);
end

%% ------------------------------------------------------------
%  20) EXPORT geometry for ANSYS
%  ------------------------------------------------------------
qsm_selected = QSM_simple;
n_opt     = length('QSM_simple');
geom_orig = myfun.result_ansys(QSM_simple, n_opt);

for i = 1:n_opt
    geom_table = geom_orig{i};                            % table of one model
    file_name  = sprintf('%s%d.txt', export_prefix, i);   % e.g. geom_IND07_v3_1.txt
    writematrix(geom_table, file_name, 'Delimiter', '\t');
    fprintf('Exported: %s\n', file_name);
end