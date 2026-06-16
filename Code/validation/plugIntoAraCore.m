% Plug a species-specific version of the pcm into AraCore
clearvars -except gurobiAvailable projDir; clc;

% specify the species that should be used to plug into AraCore. Leave empty
%   to use the union PCM. Specify the version of the pcm that should be
%   used
useSpecies = 'Arabidopsis_thaliana';%'Nicotiana_tabacum'; % 
importVersion = 2;

useUnion = isempty(useSpecies);

%% load Aracore model
% read the .mat because the grRules and KEGG associations are
% retained better than if formatting back to .xml
cyt = readCbModel( ...
    [projDir 'Data/analysis/AraCore_v2_1.wKEGG.mat']);

%% read pcm and plug into aracore
if useUnion
    modelPath = strcat(projDir, 'Data/pcm/');
    load(strcat(modelPath, 'pcm.v', num2str(importVersion), '.mat'), 'model');
else
    modelPath = strcat(projDir, 'Data/pcm/species/');
    load(strcat(modelPath, useSpecies, '.pcm.v', num2str(importVersion), '.mat'), 'model');
end
chl = model;

% sequence of inputs for the model used int the paper: 0, 0, 0, 1, 1, 1, 1
[model, missingMets] = plugAndPlay(cyt, chl);

%% postprocessing:
model = buildRxn2subSystem(model, false);
% light was not merged automatically, as apparent from missingMets
annoFieldsList = {'metCharges', 'metFormulas', 'metSEEDID', ...
    'metisinchikeyID', 'metKEGGID'};
model = mergeMets(model, 'C00205[h]', 'hnu[h]', annoFieldsList);

% Mg2+ is essential for pigment synthesis, which is missing in Aracore.
% Still add an import for Mg
model = addReaction(model, 'Im_Mg', 'reactionName', ...
    'Import Mg2+', 'reactionFormula', '--> C00305[c]');
model = addToSubsystem(model, {'Im_Mg'}, 'import');

% Aracore assumed AT4G22930, Dihydroorotase to be active in the plastid,
% but newer findings place it in the cytosol. This reaction is thus
% missing from the pcm. Add it to the cytosol and adapt transports
% remove reactions: Tr_DHO1
model = removeRxns(model, {'Tr_DHO1'});
% add DHOase_m based on previous AraCore implementation in chloroplast
model = addReaction(model, 'DHOase_m', 'reactionName', 'DHOase [m]', ...
    'reactionFormula', 'CAs[m] --> DHO[m] + H2O[m]', 'geneRule', ...
    'AT4G22930');
model = addToSubsystem(model, {'DHOase_m'}, 'Pyrimidine metabolism');
% also add CAs transporter
model = addReaction(model, 'Tr_CAs_m', 'reactionName', 'transport, CAs',...
    'reactionFormula', 'CAs[m] <=> CAs[c]');
model = addToSubsystem(model, {'Tr_CAs_m'}, 'Pyrimidine metabolism');
model = addToSubsystem(model, {'Tr_CAs_m'}, 'transport');
model = updateFromGrRules(model);
% remove Mets DHO[h]
%model = removeMetabolites(model, {'DHO[h]'}, true, 'exclusive');

% write to disk
if isfield(model, 'A')
    model = rmfield(model, 'A');
end

writeCbModel(model, 'fileName', ...
    [projDir 'Data/pcm/pluggedIntoAraCore/AraCore.' useSpecies '.pcm.xml']);
writeCbModel(model, 'fileName', ...
    [projDir 'Data/pcm/pluggedIntoAraCore/AraCore.' useSpecies '.pcm.mat']);
