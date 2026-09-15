%% Analyze models with different thresholds
clearvars -except gurobiAvailable projDir; clc;

analysisPath = [projDir 'Data/analysis/panGenomeAnalysis/' ...
    'demand_threshold_sensitivity/'];

% load the union model for reference for the missingInOrg vector
importVersion = 1;
load(strcat(projDir, 'Data/pcm/pcm.v', num2str(importVersion), '.mat'), ...
    'model');
unionModel = model;

% read media definitions
media = readtable([projDir 'Data/analysis/simulationMedia.tsv'], ...
    'FileType', 'text', 'Delimiter', '\t');
media.hetero = media.hetero_str;
media = removevars(media, ["hetero_suc", "hetero_str"]);

% calculate demand production rates for union model
[UnionPhotoProd, unionHeteroProd] = calcDemandRate(unionModel, media);

% gather list of all orgs
fileList = dir(fullfile([projDir 'Data/sequences/'], '*.fasta'));
orgs = cell(length(fileList), 1);
for i = 1:length(fileList)
    % Remove the suffix to get the base name
    [~, baseName, ~] = fileparts(fileList(i).name);
    orgs{i} = baseName;
end

% initialize tested thresholds
thresholds = 0.0:0.1:1.0;

% initialize results data structure
nModels = length(orgs) * length(thresholds);
data.threshold = zeros(nModels, 1);
data.species = cell(nModels, 1);
data.nRxns = zeros(nModels, 1);
data.nRxnsWithSupport = zeros(nModels, 1);
data.nRxnsWithoutSupport = zeros(nModels, 1);
data.nRxnsDiscarded = zeros(nModels, 1);
data.photoProdRatio = zeros(nModels, 1);
data.heteroProdRatio = zeros(nModels, 1);

% iterate over thresholds and species
modelIx = 0;
for t = 0.0:0.1:1.0
    disp(['Working on ' num2str(t)])
    rxnPAV = [];
    orgsWorked = zeros(length(orgs));
    for orgIx = length(orgs):-1:1
        org = orgs{orgIx};
        modelIx = modelIx + 1;

        data.threshold(modelIx) = t;
        data.species{modelIx} = org;

        % load missingInOrg and resulting carved model
        currModelPath = [analysisPath '04_carved_models/' org '.' num2str(t) 'pcm.v1.mat'];
        if ~isfile(currModelPath)
            % this case happened 4 times, in those, no solution was found
            % by the solver. It will appear as all-0 in the data
            continue
        end
        orgsWorked(orgIx) = 1;
        load(currModelPath)
        load([analysisPath '01_genetic_evidence_vectors/' org '.mat'])
        rxnsRetained = ismember(unionModel.rxns, model.rxns);
        [photoProd, heteroProd] = calcDemandRate(model, media);

        % record interesting numbers
        % number of reactions
        data.nRxns(modelIx) = length(model.rxns);
        % number of reactions with organism-specific genes (1 in missingInOrg)
        data.nRxnsWithSupport(modelIx) = sum(rxnsRetained & ~missingInOrg);
        % number of reactions retained for biomass-generation or flux consistency
        data.nRxnsWithoutSupport(modelIx) = sum(rxnsRetained & missingInOrg);
        % number of reactions with evidence discarded for flux consistency reasons
        data.nRxnsDiscarded(modelIx) = sum(~rxnsRetained & ~missingInOrg);
        % actual demand production (ratio)
        data.photoProdRatio(modelIx) = photoProd / UnionPhotoProd;
        data.heteroProdRatio(modelIx) = heteroProd / unionHeteroProd;

        % PAV 
        rxnPAV = [rxnPAV; transpose(rxnsRetained)];
    end
    
    % output resulting rxnPav
    resTable = array2table(rxnPAV, 'VariableNames', unionModel.rxns, ...
        'RowNames', orgs(logical(orgsWorked)));
    outPath = [analysisPath 'rxnPav_' num2str(t) '.csv'];
    writetable(resTable, outPath, "WriteRowNames",true);
end

% write data to disk
table = struct2table(data);
writetable(table, [analysisPath 'sensitivityData.csv']);


%% Functions
function [photoProd, heteroProd] = calcDemandRate(model, media)
    oriUb = model.ub;

    % Phototrophic conditions
    for rxnIx = 1:length(media.rxn)
        currRxnIx = strcmp(media.rxn{rxnIx}, model.rxns);
        if any(currRxnIx)
            model.ub(currRxnIx) = media.photo(rxnIx);
        end
    end
    photoProd = solveLP(model).f;    
    
    % Heterotrophic conditions
    for rxnIx = 1:length(media.rxn)
        currRxnIx = strcmp(media.rxn{rxnIx}, model.rxns);
        if any(currRxnIx)
            model.ub(currRxnIx) = media.hetero(rxnIx);
        end
    end
    heteroProd = solveLP(model).f;
    model.ub = oriUb;
end