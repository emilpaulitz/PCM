% This script extracts the kcat values from the enzyme-constrained models
clearvars -except gurobiAvailable projDir; clc;

dataPath = [projDir 'Code/gecko/panAraEcPcm/data/'];
inModelPath = [projDir 'Data/analysis/accSpecPCM/04_models/'];
outModelPath = [projDir 'Code/gecko/panAraEcPcm/models/'];

fileList = dir(fullfile(inModelPath, '*.mat'));
media = readtable([projDir 'Data/analysis/simulationMedia.tsv'], ...
    'FileType', 'text', 'Delimiter', '\t');

% laodEcModel requires a model Adapter to be set and takes the gecko dir
% path from it. That one is the same across all accession-adapters
adapterLocation = fullfile(projDir, ['Code/gecko/panAraEcPcm/' ...
    'adapters/Col_0/panAraEcPcmAdapter.m']);
ModelAdapterManager.setDefault(adapterLocation);

outputTables = {};
for accIx = 1:length(fileList)
    fileName = fileList(accIx).name;
    acc = erase(fileName, '.mat');
    accUnderscore = strrep(acc, '-', '_');
    disp(['Starting accession ' acc '!'])

    ecModel = loadEcModel(['ath.' accUnderscore '.ecModel.yml']);

    rxnsData = ecModel.ec.rxns;
    kcatData = ecModel.ec.kcat;

    % Create a table for current accession
    accTable = table(rxnsData, kcatData, 'VariableNames', ...
        {[acc '_rxns'], [acc '_kcat']});

    % Add table to the list for concatenation
    outputTables{end+1} = accTable; 
end

% Write each table in outputTables to its own CSV file beause row numbers
% do not match
for i = 1:length(outputTables)
    currentTable = outputTables{i};
    outputFileName = fullfile(projDir, ...
        sprintf('Data/analysis/accSpecPCM/05_comparison_data/kcats/%s_kcats.csv', ...
                erase(fileList(i).name, '.mat')));
    
    writetable(currentTable, outputFileName);
end