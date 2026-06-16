clearvars -except gurobiAvailable projDir; clc;

% load pcm
dataPath = [projDir 'Data/analysis/cyanoCompare/'];
modelPath = [projDir 'Data/pcm/'];
load(strcat(modelPath, 'pcm.v1.mat'), 'model');
pcm = model;

% load cyanobacterium, curated with modelCuration.m
cya = readCbModel([dataPath 'iSynCJ816_curated.xml']);
% gene 205 is a duplication of 204 and causes the rule to be faulty
cya.rules(strcmp(cya.rxns, 'GLUt4pp')) = {'(( x(203) & x(204) & x(206) ) | x(584) )'};
cya.genes{205} = 'sll1103and';
cya = removeGenesFromModel(cya, 'sll1103and');

% add KEGG IDs 
cya_mets_kegg = readtable([dataPath 'cya_to_kegg.csv'], 'FileType', 'text', 'Delimiter', ',');
cya.metKEGGID = cell(size(cya.mets));
for i = 1:length(cya.mets)
    matchIdx = find(strcmp(cya_mets_kegg.cyano_met_id, cya.mets{i}), 1);
    
    if ~isempty(matchIdx)
        cya.metKEGGID{i} = cya_mets_kegg.kegg_id{matchIdx};
    else
        cya.metKEGGID{i} = '';
    end
end

%% load AraCode with metKEGGID
cyt = readCbModel([projDir 'Data/analysis/AraCore_v2_1.wKEGG.mat']);
% generate subSystemNames and remove space after sep comma
newSubsystemNames = cell.empty;
for i = 1:length(cyt.subSystems)
    newSubsystemNames = unique([newSubsystemNames, ...
        strsplit(cyt.subSystems{i}, ', ')]);
    cyt.subSystems{i} = strjoin(strsplit(cyt.subSystems{i}, ', '), ',');
end
cyt.subSystemNames = newSubsystemNames;
% generate rxn2subsystems
cyt.rxn2subSystem = zeros(length(cyt.rxns), length(newSubsystemNames));
for i = 1:length(cyt.rxns)
    for j = 1:length(cyt.subSystemNames)
        if any(strcmp(cyt.subSystems{i}, cyt.subSystemNames{j}))
            cyt.rxn2subSystem(i, j) = 1;
        end
    end
end

%% Preprocess the cyanobacterium model
% change the name of the cyanobacterium cytosol to stroma, and its
% compartment symbol to h, to avoid confusion
cya.comps{1} = 'h';
cya.compNames{1} = 'Chloroplast';
cya.mets = cellfun(@(met) strrep(met, 'c[c]', 'h[h]'), cya.mets, ...
    'UniformOutput', false);
% also rename the cyanobacterium extracellular to cytosol, to be able to
% match the whole-cell cytoplasm
cya.comps{2} = 'c';
cya.compNames{2} = 'Cytosol';
cya.mets = cellfun(@(met) strrep(met, 'e[e]', 'c[c]'), cya.mets, ...
    'UniformOutput', false);

% change format of subsystems to match the pcm's
newSubSystem = cell(size(cya.subSystems));
newSubsystemNames = cell.empty;
for i = 1:length(cya.subSystems)
    newSubSystem{i} = strjoin(cya.subSystems{i}, ',');
    currSubs = cya.subSystems{i};
    newSubsystemNames = unique([newSubsystemNames; cya.subSystems{i}]);
end    
cya.subSystemNames = newSubsystemNames;
% generate rxn2subsystems
cya.rxn2subSystem = zeros(length(cya.rxns), length(newSubsystemNames));
for i = 1:length(cya.rxns)
    for j = 1:length(cya.subSystemNames)
        if any(strcmp(cya.subSystems{i}, cya.subSystemNames{j}))
            cya.rxn2subSystem(i, j) = 1;
        end
    end
end
cya.subSystems = newSubSystem;

% then put all exchange/export/import rxns in their respective subsystems
for i = 1:length(cya.rxns)
    if startsWith(cya.rxns{i}, 'EX_')
        if cya.ub(i) == 0
            cya = addSubs(cya, i, 'export');
        else
            if cya.lb(i) == 0
                cya = addSubs(cya, i, 'import');
            else
                cya = addSubs(cya, i, 'exchange');
            end
        end
    end
end

% add grRules field
cya = creategrRulesField(cya);

% rename/add some fields to match pcm naming convention
cya.rxnEC = cya.rxnECNumbers;
cya.rxnNotes = repmat({''},size(cya.rxns));
cya.rxnKEGGID = repmat({''},size(cya.rxns));

%% Plug it in
% series of inputs: 0, 1, [3x] (enter), z, [4x] 0, z, 2
[model, missingMets] = plugAndPlay(cyt, cya);

%% verify and curate result
%verifyModel(model) % this outputs subsystem formatting inconsistencies
%missingMets % this outputs light and starch

% merge photons; They were in different compartments
metAnnoFields = {'metKEGGID', 'metSEEDID', 'metisinchikeyID'};
model = mergeMets(model, 'photon_c[c]', 'hnu[h]', metAnnoFields);

% replace starch-2 by double the amount of glycogen in the biomass reaction
bmIx = find(model.S(strcmp(model.mets, 'starch2[h]'), :));
model.S(strcmp(model.mets, 'glycogen_h[h]'), bmIx) = ...
    2 * model.S(strcmp(model.mets, 'starch2[h]'), bmIx);
model.S(strcmp(model.mets, 'starch2[h]'), bmIx) = 0;

% merge sodium ion metabolites; They are already matched in the plastid
model = mergeMets(model, 'na1_c[c]', 'Na[c]', metAnnoFields);

%% to be able to compare with the pcm version, also add the mitochondrial
% dihydroorotase. This should have been mitochondrial in AraCore, actually
model = removeRxns(model, {'Tr_DHO1'});
% add DHOase_m based on previous AraCore implementation in chloroplast
model = addReaction(model, 'DHOase_m', 'reactionName', 'DHOase [m]', ...
    'reactionFormula', 'CAs[m] --> DHO[m] + H2O[m]', 'geneRule', ...
    'AT4G22930');
% also add CAs transporters
model = addReaction(model, 'Tr_CAs_m', 'reactionName', 'transport, CAs [m]',...
    'reactionFormula', 'CAs[m] <=> CAs[c]');
model = addToSubsystem(model, {'Tr_CAs_m'}, 'transport');
model = addReaction(model, 'Tr_CAs_h', 'reactionName', 'transport, CAs [h]',...
    'reactionFormula', 'cbasp_h[h] <=> CAs[c]');
model = addToSubsystem(model, {'Tr_CAs_h'}, 'transport');
model = updateFromGrRules(model);

%% output merged model
% output as mat
if isfield(model, 'A')
    pcm = rmfield(model, 'A');
end
writeCbModel(model, 'fileName', strcat(dataPath, 'aracore.cyano.mat'));

% output as xml
if ~isfield(model, 'rxnECNumbers') && isfield(model, 'rxnEC')
    model.rxnECNumbers = strrep(model.rxnEC, ',', ';');
end

% Change format of the subsystems array
% This will give warnings with writeCbModel, but will 
% map subsystems correctly in the .xml and for plotting purposes
% subsNames and rxn2subs are not present in Aracore, so skip this for now.
for i = 1:length(model.rxns)
    %s = model.subSystemNames(logical(model.rxn2subSystem(i, :)));
    %model.subSystems{i} = s;
end

% write model
writeCbModel(model, 'format','sbml', 'filename', ...
    strcat(dataPath, 'aracore.cyano'));

%% functions
% first checks whether the subsystem exists, and then calls the respective
% function
function model = addSubs(model, rxnIx, subs)
    if any(strcmp(subs, model.subSystemNames))
        model = addToSubsystem(model, model.rxns(rxnIx), subs);
    else
        model = addSubSystem(model, model.rxns(rxnIx), subs);
    end
end
