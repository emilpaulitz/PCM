% This script was run on a machine with Gurobi available
clearvars -except gurobiAvailable projDir; clc;

importVersion = 2;

outputDir = [projDir 'Data/analysis/accSpecPCM/05_comparison_data/'];
modelsDir = [projDir 'Data/analysis/accSpecPCM/04_models/'];

% load union PCM model and obtain set of all rxns
modelPath = strcat(projDir, 'Data/pcm/');
load(strcat(modelPath, 'pcm.v', num2str(importVersion), '.mat'), 'model');
pcm = model;
allRxns = sort(model.rxns);

% Get the list of model files from modelsDir
modelFiles = dir(strcat(modelsDir, '*.mat'));

% Load all accession-specific models into a map
models = containers.Map();
for file = modelFiles'
    baseName = file.name(1:end-length('.mat')); 

    % Add model to the map with the base name as key
    modelData = load(strcat(file.folder, '/', file.name), 'model');
    models(baseName) = modelData.model; 
end

% Initialize output struct
nModels = length(models);
data.model_name = keys(models)'; % Extract the names of models from the map
data.nRxns = zeros(nModels, 1);
data.nMets = zeros(nModels, 1);

% Gather information from each model
for modelIx = 1:length(data.model_name)
    currModel = models(data.model_name{modelIx});
    data.nRxns(modelIx) = length(currModel.rxns);
    data.nMets(modelIx) = length(currModel.mets);
end

% Output to file
csvFilePath = fullfile(outputDir, 'modelData.csv');
writetable(struct2table(data), csvFilePath);

% rxn PAV - Create a presence/absence matrix
columnNames = allRxns;
resArray = [];
accs = keys(models);
for modelName = accs
    currModel = models(modelName{1});
    orgRxns = sort(currModel.rxns);
    resArray = [resArray; transpose(ismember(allRxns, orgRxns))];
end

% Output resulting rxnPAV
resTable = array2table(resArray, 'VariableNames', columnNames, 'RowNames', accs);
outPath = [outputDir 'rxnPav.csv'];
writetable(resTable, outPath, "WriteRowNames", true);
