% Gather information about reaction presence absence in the pcm versions
clearvars -except gurobiAvailable projDir; clc;

modelPath = strcat(projDir, 'Data/pcm/');
importVersion = 1;

% load the union pcm for the set of all rxns
load([modelPath 'pcm.v' num2str(importVersion) '.mat'], 'model');
allRxns = sort(model.rxns);

% initialize result table
columnNames = allRxns;
resArray = [];

% iterate organisms and add their columns
fileList = dir(fullfile([modelPath 'species/'], ...
    ['*.pcm.v' num2str(importVersion) '.mat']));
orgs = cell(length(fileList), 1);
for fileIx = 1:length(fileList)
    currFile = fileList(fileIx).name;
    [~, baseName, ~] = fileparts(currFile);
    orgs{fileIx} = baseName;
    
    load([modelPath 'species/' currFile], 'model');
    orgRxns = sort(model.rxns);
    resArray = [resArray; transpose(ismember(allRxns, orgRxns))];

    % also write the models' genes for further analysis in Python
    if contains(baseName, 'Marchanta_polymorpha')
        baseName = ['Marchantia_polymorpha.pcm.v' num2str(importVersion)];
    end
    fileID = fopen(['Data/analysis/panGenomeAnalysis/modelGenes/' baseName '.genes.txt'], 'w');
    for i = 1:length(model.genes)
        fprintf(fileID, '%s\n', model.genes{i}); 
    end
    fclose(fileID);
end

% output resulting rxnPav
resTable = array2table(resArray, 'VariableNames', columnNames, 'RowNames', orgs);
outPath = [projDir 'Data/analysis/panGenomeAnalysis/rxnPav.csv'];
writetable(resTable, outPath, "WriteRowNames",true);
