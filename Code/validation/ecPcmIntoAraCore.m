% Plug a species-specific enzyme-constrained version of the pcm into the
% enzyme-constrained AraCore
clearvars -except gurobiAvailable projDir; clc;

modelPath = strcat(projDir, 'Data/pcm/ecPCM/');
inputPcm = 'Ath.tuned.xml';
outputFname = 'AraCore.ecpcm.tuned';

%% load Aracore model
cyt = readCbModel([projDir 'Data/analysis/ecAraCore.wKEGG.mat']);

%% read ecPcm
chl = readCbModel([modelPath inputPcm]);
% all original metabolites got an extra compartment tag after read-in
% match to compartment tag with anything inside except a square bracket
pattern = '(\[[^\]]+\])\1'; 
for i = 1:length(chl.mets)
    chl.mets{i} = regexprep(chl.mets{i}, pattern, '$1');
end
% these fields got lost somehow
chl = buildRxn2subSystem(chl, false);
chl = creategrRulesField(chl);

%% merge
% sequence: 0 (no backup) 0 (no other symbs) 0 (not remove transp)
% 0 (not add KEGGs) 1 (match anyway) 1 (match anyway) 1 (use rxn without
% additional proton)
[model, missingMets] = plugAndPlay(cyt, chl, true);

%% postprocessing:
% rxnECNumbers field was not properly retained from ecModel. We do not need
% it currently, however
if isfield(model, 'rxnECNumbers')
    model = rmfield(model, 'rxnECNumbers');
end

% light was not merged automatically, as apparent from missingMets
annoFieldsList = {'metCharges', 'metFormulas', ...
    'metisinchikeyID', 'metKEGGID'};
model = mergeMets(model, 'C00205[h]', 'hnu[h]', annoFieldsList);

% Mg2+ is essential for pigment synthesis, which is missing in Aracore.
% Still add an import for Mg
model = addReaction(model, 'Im_Mg', 'reactionName', ...
    'Import Mg2+', 'reactionFormula', '--> C00305[c]', ...
    'subSystem', 'import');

% Aracore assumed AT4G22930, Dihydroorotase to be active in the plastid,
% but newer findings place it in the cytosol. This reaction is thus
% missing from the pcm. Add it to the cytosol and adapt transports
% remove reactions: Tr_DHO1
model = removeRxns(model, {'Tr_DHO1'});
% add DHOase_m based on previous AraCore implementation in chloroplast
model = addReaction(model, 'DHOase_c', 'reactionName', 'DHOase [c]', ...
    'reactionFormula', 'CAs[c] --> DHO[c] + H2O[c]', 'geneRule', ...
    'AT4G22930');
model = updateFromGrRules(model);
% remove Mets DHO[h]
%model = removeMetabolites(model, {'DHO[h]'}, true, 'exclusive');

% write to disk
if isfield(model, 'A')
    model = rmfield(model, 'A');
end

% this speeds up file writing because the rxn2subsystems field does not
% have to be generated multiple times
model = buildRxn2subSystem(model);

writeCbModel(model, 'fileName', [projDir 'Data/pcm/pluggedIntoAraCore/' outputFname '.xml']);
writeCbModel(model, 'fileName', [projDir 'Data/pcm/pluggedIntoAraCore/' outputFname '.mat']);