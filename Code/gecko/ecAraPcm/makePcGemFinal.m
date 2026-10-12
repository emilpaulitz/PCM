% Make the Arabidopsis version of the PCM enzyme-constrained
clearvars -except gurobiAvailable projDir; clc;

modelPath = [projDir 'Data/pcm/'];
load(strcat(modelPath, 'species/Arabidopsis_thaliana.pcm.v3.mat'), 'model');

pcm = model;

% format subsystems and rm A field to write a functioning sbml
for i = 1:length(pcm.rxns)
    s = pcm.subSystemNames(logical(pcm.rxn2subSystem(i, :)));
    pcm.subSystems{i} = s;
end
if isfield(pcm, 'A')
    pcm = rmfield(pcm, 'A');
end
writeCbModel(pcm, 'fileName', ...
    [projDir, 'Code/gecko/ecAraPcm/models/ath.pcm.v3.xml']);

%% check RAVEN and GECKO
checkInstallation;
GECKOInstaller.install;

%% Set default model adapter
% This defines the model version that is loaded!
adapterLocation = fullfile(projDir, 'Code/gecko/ecAraPcm', ...
    'ecAraPcmAdapter.m');
ModelAdapterManager.setDefault(adapterLocation);

%% load the conventional pcm
model = loadConventionalGEM();

%% create draft full ecGEM
% The data was modified to include unreviewed proteins and exclude
% duplications

dataPath = ['/home/emil/Desktop/PhD-Synch/pan_chloroplast/' ...
    'pan_chl_model/Code/gecko/ecAraPcm/data/'];
fname = [dataPath 'uniprot.tsv'];

% create draft ecModel
model.metNames(strcmp(model.mets, 'C00111[h]')) = {'glycerone phosphate'};
[ecModel, noUniprot] = makeEcModel(model);

%% gather eccodes
% correct some issues in formatting of the EC numbers
ecModel.eccodes = cellfun(@(s) regexprep(strrep(strrep(s, ' ', ''), ...
    ',', ';'), '^;', ''), ecModel.eccodes, 'UniformOutput', false);
ecModel = getECfromGEM(ecModel);
noEC = cellfun(@isempty, ecModel.ec.eccodes);
ecModel = getECfromDatabase(ecModel, noEC);

%% query BRENDA
kcatList_fuzzy = fuzzyKcatMatching(ecModel);

%% DLKcat
[ecModel, noSMILES] = findMetSmiles(ecModel);

% have to manually delete the file. The file gets locked somehow after
% running once.
writeDLKcatInput(ecModel);

% run DLKcat
runDLKcat();

% read output
kcatList_DLKcat = readDLKcatOutput(ecModel);

%% Merge BRENDA and DLKcat outputs
kcatList_merged = mergeDLKcatAndFuzzyKcats(kcatList_DLKcat, ...
    kcatList_fuzzy);

% and assign to ecModel
ecModel = selectKcatValue(ecModel, kcatList_merged);

% use average kcats for isozymes (otherwise, missing kcats will allow very
% high rates)
ecModel = getKcatAcrossIsozymes(ecModel);

% use average kcat and MW for missing enzymes
[ecModel, rxnsMissingGPR, standardMW, standardKcat] = ...
    getStandardKcat(ecModel);

% apply to ecModel
ecModel = applyKcatConstraints(ecModel);


%% set ptot
ecModel = setProtPoolSize(ecModel);

%% Tune too low kcats until Rubisco would be the next protein to be tuned
% Before, used the tunedKcats output for reverting the changes. However,
% that table is not sorted by tuning iterations but by name of the reaction
% since some enzymes might be tuned multiple times. Therefore, implement
% this procedure by changing the code in sensitivityTuning.m (Look for
% "% Custom code start" in sensitivityTuningRubisco.m)
ecModel = setProtPoolSize(ecModel);
[ecModel, tunedKcats] = sensitivityTuningRubisco(ecModel);
%struct2table(tunedKcats)

% apply custom kcats
[ecModel, rxnUpdated, notMatch] = applyCustomKcats(ecModel);

% apply new Kcats
ecModel = applyKcatConstraints(ecModel);
ecModel = setProtPoolSize(ecModel);

%% Write to disk
% Writing to and reading from yaml not only introduces errors but shuffles
% Miriams! Do not do it!
%saveEcModel(ecModel, 'ath.ecModel.yml');
%ecModel_yml = loadEcModel('ath.ecModel.yml');

ecModel.id = 'Ath ecPCM';
exportModel(ecModel, [modelPath 'ecPCM/Ath.tuned.xml']);
