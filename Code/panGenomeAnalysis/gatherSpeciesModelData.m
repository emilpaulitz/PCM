%% Analyze the size of resulting species-specific models
clearvars -except gurobiAvailable projDir; clc;

% Get list of files 
folderPath = [projDir 'Data/pcm/species/'];
matFiles = dir(fullfile(folderPath, '*v2.mat'));

% Load models
modelMap = containers.Map;
disp('Loading models...')
for k = 1:length(matFiles)
    matFilePath = fullfile(folderPath, matFiles(k).name);
    
    loadedModel = load(matFilePath, 'model');
    
    % Check if the 'model' variable was correctly loaded
    if isfield(loadedModel, 'model')
        
        % Trim '.pcm.v2.mat' from the filename to use as the key in the map
        [~, fileName, ~] = fileparts(matFiles(k).name);
        key = fileName(1:end-length('.pcm.v2'));
        
        % Add the model to the map with the species name as key
        modelMap(key) = loadedModel.model;
    else
        disp(['Warning: "model" variable not found in file: ', matFiles(k).name]);
    end
end

% Get size of the models
species = keys(modelMap);
keyNames = cell(length(species), 1);
fieldLengths = zeros(length(species), 1);
for i = 1:length(species)
    key = species{i};
    model = modelMap(key);
    
    % Check if the specified field exists in the model
    if isfield(model, 'rxns')
        keyNames{i} = key;
        fieldLengths(i) = length(model.rxns);
    else
        disp(['Warning: Field "', fieldName, '" not found in model associated with key: ', key]);
    end
end

% Create a table with key names and field lengths
resultTable = table(keyNames, fieldLengths, ...
    'VariableNames', {'Species', 'nRxns'});
sortrows(resultTable, 'nRxns')