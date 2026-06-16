% Make the Arabidopsis version of the PCM enzyme-constrained
clearvars -except gurobiAvailable projDir; clc;

modelPath = [projDir 'Data/pcm/'];
load(strcat(modelPath, 'species/Arabidopsis_thaliana.pcm.v2.mat'), 'model');
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
    [projDir, 'Code/gecko/ecAraPcm/models/ath.pcm.v2.xml']);

%% check RAVEN and GECKO
checkInstallation;
GECKOInstaller.install;

%% Set default model adapter
adapterLocation = fullfile(projDir, 'Code/gecko/ecAraPcm', ...
    'ecAraPcmAdapter.m');
ModelAdapterManager.setDefault(adapterLocation);

%% load the conventional pcm
model = loadConventionalGEM();

%% create draft full ecGEM
% The data was modified to include unreviewed proteins and exclude
% duplications. See makePcGem.m

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

%% Tune too low kcats until Rubisco is the most used protein
ecModel = setProtPoolSize(ecModel);
[ecModel, tunedKcats] = sensitivityTuning(ecModel);
%struct2table(tunedKcats)

% Revert back all kcats that came after rubisco, including rubisco
rbcIx = find(startsWith(tunedKcats.rxns, 'RBC_h') | ...
    startsWith(tunedKcats.rxns, 'RBO_h'), 1 );
for i = rbcIx:length(tunedKcats.rxns)
    currRxnIx = strcmp(ecModel.ec.rxns, tunedKcats.rxns{i});
    ecModel.ec.kcat(currRxnIx) = tunedKcats.oldKcat(i);
end

% apply custom kcats
[ecModel, rxnUpdated, notMatch] = applyCustomKcats(ecModel);

% apply new Kcats
ecModel = applyKcatConstraints(ecModel);
ecModel = setProtPoolSize(ecModel);

%% Write to disk
saveEcModel(ecModel, 'ath.ecModel.yml');
ecModel = loadEcModel('ath.ecModel.yml');

ecModel.id = 'Ath ecPCM';
exportModel(ecModel, [modelPath 'ecPCM/Ath.tuned.xml']);
