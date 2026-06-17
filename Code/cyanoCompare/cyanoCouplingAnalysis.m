% Performs flux coupling analysis with the F2C2 package
% Please install F2C2 by downloading it and extracting into 
% [projDir 'Resources/F2C2 v0.95b/']

clearvars -except gurobiAvailable projDir; clc;
addpath([projDir 'Resources/F2C2 v0.95b/code'])

modelsPath = [projDir 'Data/analysis/cyanoCompare/coupling_models/'];
outPath = [projDir 'Data/analysis/cyanoCompare/coupling_results/'];

% Read models
% Read AraCore
ara = readCbModel([projDir 'Data/analysis/AraCore_v2_1.wKEGG.mat']);
% add a night AraCore that uses starch2
nightAra = ara;
nightAra.ub(strcmp(nightAra.rxns, 'Im_hnu')) = 0;
nightAra.ub(strcmp(nightAra.rxns, 'Im_CO2')) = 0;
nightAra.lb(strcmp(nightAra.rxns, 'Ex_O2')) = -1000;
nightAra = addReaction(nightAra, 'Im_starch', 'reactionFormula', ...
                        '--> starch2[h]', 'lowerBound', 0, 'upperBound', 8);

% Read APC
dayApc = readCbModel([modelsPath 'day_apc_trin.xml']);
nightApc = readCbModel([modelsPath 'night_apc_trin.xml']);
% Read ACY
dayAcy = readCbModel([modelsPath 'day_acy_trin.xml']);
nightAcy = readCbModel([modelsPath 'night_acy_trin.xml']);

%% Apply F2C2
models = {ara, nightAra, dayApc, nightApc, dayAcy, nightAcy};
modelNames = {'AraCore', 'nightAraCore', 'dayApc_trin', 'nightApc_trin',...
              'dayAcy_trin', 'nightAcy_trin'};
for modelIx = 1:length(models)
    disp(['Working on ' modelNames{modelIx}]);

    conformModel = models{modelIx};

    % flip reversed reactions so that all reactions are either reversible 
    % or irreversible from left to right
    for i=1:length(conformModel.rxns)
        if  conformModel.lb(i) < 0 && conformModel.ub(i) <= 0
            conformModel = flipReversed(conformModel,conformModel.rxns(i));
        end
    end

    if any(conformModel.ub == 0)
        conformModel = removeReactions(conformModel, ...
            conformModel.rxns(conformModel.ub == 0));
    end

    if any(conformModel.ub < 0)
       disp('error found incorrect rxn');
       continue 
    end
    
    conformModel = makeModelComform(conformModel);
    
    % all other solvers are not usable
    [fctable, blocked] = F2C2('glpk', conformModel);

    % Write out data
    rxns = conformModel.Reactions(~blocked);
    cellData = [rxns, num2cell(fctable)];
    allColNames = [{' '}; rxns];
    T = cell2table(cellData, 'VariableNames', allColNames);
    writetable(T, [outPath modelNames{modelIx} '.csv']);
end

% Analyze data in Python
%% FUNCTIONS
function model = makeModelComform(inputModel)
    model = struct();
    model.stoichiometricMatrix = full(inputModel.S);
    model.reversibilityVector = double(inputModel.lb < 0);
    model.Reactions = inputModel.rxns;
    model.Metabolites = inputModel.mets;
end