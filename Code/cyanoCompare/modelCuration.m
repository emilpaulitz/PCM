% This script and the ensuing model "iSynCJ816_curated.xml was kindly 
% shared by Mauricio Alexander de Moura Ferreira!
%
% Manual curation on the iSynCJ816 model provided by Joshi et al., Algal 
% Research, 2017. The modifications performed here were based on the
% changes made by Andrews et al., Algal Research, 2024 in preparation for
% the sMOMENT model that was built afterwards.
%
%%
% initCobraToolbox(false);
cd('~/Posdoutorado/2_Synechocystis_GEM/Code/');
clc;clear

%% Load the model
% modelBiGG = readCbModel('/home/mferreira/Posdoutorado/2_Synechocystis_GEM/Models/iSynCJ816_BiGG.xml');
model = readCbModel('/home/mferreira/Posdoutorado/2_Synechocystis_GEM/Models/iSynCJ816_mmc12.xml');

%% Convert to RAVEN format since RAVEN functions are slighly better for handling changes to models than COBRA functions
modelRAVEN = ravenCobraWrapper(model);
modelRAVEN.rxnReferences = model.rxnReferences;

%% Fixing RBPC
% model.S(369, 827) = 1;

%% Add reactions to model
rxnsToAdd = struct();

% Adding cytosolic HCO3 equilibration reaction
rxnsToAdd.rxns{1,1} = 'HCO3E2';
rxnsToAdd.rxnNames{1,1} = 'HCO3 equilibration reaction (cytosol)';
rxnsToAdd.equations{1,1} = 'hco3_c + h_c => h2o_c + co2_c';
rxnsToAdd.subSystems{1,1} = 'Carbon fixation';
rxnsToAdd.grRules{1,1} = '';

% Add isoprene production and exchange reaction
metsToAdd = struct();
metsToAdd.mets = {'isoprene_c', 'isoprene_e', 'isoprene_p'}';
metsToAdd.metNames = {'Isoprene', 'Isoprene', 'Isoprene'}';
metsToAdd.compartments = {'c', 'e', 'p'}';

rxnsToAdd.rxns{2,1} = 'ISPS';
rxnsToAdd.rxnNames{2,1} = 'Isoprene synthase';
rxnsToAdd.equations{2,1} = 'dmpp_c => isoprene_c';
rxnsToAdd.subSystems{2,1} = 'Terpenoid backbone biosynthesis';
rxnsToAdd.grRules{2,1} = 'IspS';

rxnsToAdd.rxns{3,1} = 'ISPtpp';
rxnsToAdd.rxnNames{3,1} = 'Isoprene transport (periplasm)';
rxnsToAdd.equations{3,1} = 'isoprene_p <=> isoprene_c';
rxnsToAdd.subSystems{3,1} = 'Transport';
rxnsToAdd.grRules{3,1} = '';

rxnsToAdd.rxns{4,1} = 'ISPtex';
rxnsToAdd.rxnNames{4,1} = 'Isoprene transport (extracellular to periplasm)';
rxnsToAdd.equations{4,1} = 'isoprene_e <=> isoprene_p';
rxnsToAdd.subSystems{4,1} = 'Transport';
rxnsToAdd.grRules{4,1} = '';

rxnsToAdd.rxns{5,1} = 'EX_isoprene(e)';
rxnsToAdd.rxnNames{5,1} = 'Isoprene exchange';
rxnsToAdd.equations{5,1} = 'isoprene_e <=> ';
rxnsToAdd.subSystems{5,1} = 'Terpenoid backbone biosynthesis';
rxnsToAdd.grRules{5,1} = '';

genesToAdd = struct();
genesToAdd.genes = 'IspS';
modelRAVEN = addGenesRaven(modelRAVEN, genesToAdd);

% Add biomass exchange reaction
metsToAdd.mets{4,1} = 'biomass_c';
metsToAdd.metNames{4,1} = 'Biomass';
metsToAdd.compartments{4,1} = 'c';

modelRAVEN = addMets(modelRAVEN, metsToAdd);

modelRAVEN.S(933, 74) = 1; % add biomass pseudometabolite to autotrophic biomass reaction
modelRAVEN.S(933, 75) = 1; % add biomass pseudometabolite to mixotrophic biomass reaction
modelRAVEN.S(933, 76) = 1; % add biomass pseudometabolite to heterotrophic biomass reaction

rxnsToAdd.rxns{6,1} = 'EX_growth(c)';
rxnsToAdd.rxnNames{6,1} = 'Growth';
rxnsToAdd.equations{6,1} = 'biomass_c =>';
rxnsToAdd.subSystems{6,1} = 'Biomass';
rxnsToAdd.grRules{6,1} = '';

% Add reactions
modelRAVEN = addRxns(modelRAVEN, rxnsToAdd);

%% Rename reactions and metabolites to (somewhat) match the BiGG model
% Remove 'R_' prefix
for i = 1:length(modelRAVEN.rxns)
    modelRAVEN.rxns{i} = strrep(modelRAVEN.rxns{i}, 'R_', '');
end

% Replace (e) with _e
for i = 1:length(modelRAVEN.rxns)
    % Get the current string
    currentStr = modelRAVEN.rxns{i};
    
    % Check if the last three characters contain parentheses
    if length(currentStr) >= 3
        lastThreeChars = currentStr(end-2:end);
        
        % Check if they are of the form (X) where X is a letter
        if lastThreeChars(1) == '(' && lastThreeChars(3) == ')' && isletter(lastThreeChars(2))
            % Replace '(' with '[' and ')' with ']'
            currentStr(end-2) = '_';
            currentStr(end) = '';
        end
    end
    
    % Update the cell array with the modified string
    modelRAVEN.rxns{i} = currentStr;
end

% Rename glucose exchange
modelRAVEN.rxns{783,1} = 'EX_glc__D_e';

%% Other adjustments
% Fix subsystems
modelRAVEN.subSystems{74,1} = 'Biomass';

% Fix gene
modelRAVEN.genes{806} = 'sll1103';

% Fix subsystems
modelRAVEN.subSystems{1046,1} = char(modelRAVEN.subSystems{1046,1});
modelRAVEN.subSystems{1047,1} = char(modelRAVEN.subSystems{1047,1});
modelRAVEN.subSystems{1048,1} = char(modelRAVEN.subSystems{1048,1});
modelRAVEN.subSystems{1049,1} = char(modelRAVEN.subSystems{1049,1});
modelRAVEN.subSystems{1050,1} = char(modelRAVEN.subSystems{1050,1});
modelRAVEN.subSystems{1051,1} = char(modelRAVEN.subSystems{1051,1});

% Fix bounds
modelRAVEN.ub(modelRAVEN.ub == 999999) = 1000;
modelRAVEN.lb(modelRAVEN.lb == -999999) = -1000;

%% Sort all IDs
modelRAVEN = sortIdentifiers(modelRAVEN);

%% Export model
exportModel(modelRAVEN, '../Models/iSynCJ816_curated.xml');
% exportToExcelFormat(modelRAVEN, '../Models/iSynCJ816_curated.xlsx');
