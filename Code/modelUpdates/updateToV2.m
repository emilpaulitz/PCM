clearvars -except gurobiAvailable projDir; clc;

modelPath = strcat(projDir, 'Data/pcm/');
importVersion = 1;
load([modelPath 'pcm.v' num2str(importVersion) '.mat'], 'model');

%% General changes
disp('Working on general changes')
% Update the subsystem of import reactions
model = addToSubsystem(model, {'Im_suc', 'Im_str', 'Im_o2'}, 'import');

% G6P (beta-G6P) annotation with C00092 (general G6P) leads to problems in
% plugAndPlay because C00092 exists as its own metabolite
model.metKEGGID{strcmp(model.mets, 'G6P[h]')} = 'C01172';

% Update metComps array
for i = 1:length(model.mets)
    model.metComps(i) = find(strcmp(model.comps, extractComp(model.mets{i})));
end

% Add source for a bunch of transporters that I had curated:
model.rxnNotes(strcmp(model.rxns, 'NDH_h')) = {'10.7554/eLife.49305'};
model.rxnNotes(strcmp(model.rxns, 'Tr_xanth')) = {'10.1073/pnas.2502160122'};
model.rxnNotes(strcmp(model.rxns, 'Tr_THF')) = {'10.1111/j.1399-3054.2006.00587.x'};
model.rxnNotes(strcmp(model.rxns, 'Tr_pABA')) = {'10.1111/j.1399-3054.2006.00587.x'};

%% Changes to address stoichiometric consistency and dead-ends
disp('Addressing stoichiometric consistency and dead-ends')
rxnsToLookAt = cell(0, 1);
rxnErrorI = 0;

% These are dead-end and take place in the mitochondrion (branched-chain 
% amino acid degradation): incorrectly included in the template models
model = removeRxns(model, {'R07599', 'R07603', 'R07601'});

% remove dead-end reactions with general metabolites
model = removeRxns(model, {'R02281', 'R12571', 'R01316', 'R01313', ...
                           'R06131', 'R00645', 'R00804', 'R00034', ...
                           'R01315', 'R01314' ,'R01999', 'R02687', ...
                           'R02250', 'R02688', 'R02054', 'R02053', ...
                           'R00857', 'R08618', 'R01340', 'R01341'});

% remove metabolites that are outside any 0 reactions (these are in the
% cytosol)
model = removeMetabolites(model, {'C20694[c]', 'C12287[c]'});

% reaction from xenobiotics and carcinogenesis
model = removeRxns(model, {'R07113'});

% To be able to use FAD, let it interchange reduction equivalents with NAD
% and NADP
[model, addedRxnId, errorMsg] = addKEGGReactionNew(model, 'R09520', 'unknown', ...
    'Added for modeling purpose. Some way of activating FAD must be present', ...
    'Riboflavin metabolism');
if ~isempty(errorMsg)
    rxnErrorI = rxnErrorI + 1;
    rxnsToLookAt{rxnErrorI} = [addedRxnId ': ' errorMsg];
end
model.rxn2subSystem(strcmp(model.rxns, 'R09520'), ...
    strcmp(model.subSystemNames, 'Riboflavin metabolism')) = 1;
model.lb(strcmp(model.rxns, 'R09520')) = -1000;
model.rxnKEGGID(strcmp(model.rxns, 'R09520')) = {'R09520'};

[model, addedRxnId, errorMsg] = addKEGGReactionNew(model, 'R09748', 'unknown', ...
    'Added for modeling purpose. Some way of activating FAD must be present', ...
    'Riboflavin metabolism');
if ~isempty(errorMsg)
    rxnErrorI = rxnErrorI + 1;
    rxnsToLookAt{rxnErrorI} = [addedRxnId ': ' errorMsg];
end
model.rxn2subSystem(strcmp(model.rxns, 'R09748'), ...
    strcmp(model.subSystemNames, 'Riboflavin metabolism')) = 1;
model.lb(strcmp(model.rxns, 'R09748')) = -1000;
model.rxnKEGGID(strcmp(model.rxns, 'R09748')) = {'R09748'};

% The siroheme biosynthesis pathway's occurrence in plants is supported by 
% metacyc and occurrence of the following reactions in the chloroplast is 
% supported by Uniprot. But it is a protein, so it is out of scope for this
% model
% #Siroheme
model = removeRxns(model, {'R02864', 'R03194'});

% To connect Phosphatidate C00416[h] (its reaction has evidence in the
% chloroplast in Uniprot), add this reaction as the most direct route
% [model, addedRxnId, errorMsg] = addKEGGReactionNew(model, 'R02240', 'unknown', ...
%     ['Added to connect Phosphatidate C00416[h] since its reaction has' ...
%     'evidence in the chloroplast in Uniprot'], ...
%     'glycerolipid metabolism,glycerophospholipid metabolism');
% if ~isempty(errorMsg)
%     i = i + 1;
%     rxnsToLookAt{i} = [addedRxnId ': ' errorMsg];
% end
% model.rxn2subSystem(strcmp(model.rxns, 'R02240'), ...
%     strcmp(model.subSystemNames, 'glycerolipid metabolism')) = 1;
% model.rxn2subSystem(strcmp(model.rxns, 'R02240'), ...
%     strcmp(model.subSystemNames, 'glycerophospholipid metabolism')) = 1;
% model.rxnKEGGID(strcmp(model.rxns, 'R02240')) = {'R02240'};

% NADPH can spontaneously react to NADPHX, which is repaired by R10288. Add
% the damage reaction
model = addReaction(model, 'NADPH_dmg_h', 'reactionName', ...
    'spontaneous NADPH damaging; MNXR159750', ...
    'reactionFormula', 'C00005[h] + C00001[h] --> C04899[h]');
model.subSystems(strcmp(model.rxns, 'NADPH_dmg_h')) = {''};
% Same for NADH, repaired by R00129
model = addReaction(model, 'NADH_dmg_h', 'reactionName', ...
    'spontaneous NADH damaging; MNXR117736', ...
    'reactionFormula', 'C00004[h] + C00001[h] --> C04856[h]');
model.subSystems(strcmp(model.rxns, 'NADH_dmg_h')) = {''};

% Remove a cluster of reactions dealing with lipids, interconverting
% between general side chains and concrete metabolites. The cluster is not
% connected to the rest of the metabolism
% #phosphocholine
model = removeRxns(model, {'R04480', 'R01318', 'R02241','R01999', ...
    'R00851', 'R07064', 'R07859', 'R01317', 'R07868', 'R07864', ...
    'R07865', 'R07057', 'R08177', 'R03814', 'R01596'});

% Similar case, with dead-end Acetly-ACP
model = removeRxns(model, {'R04355'});

% Similar case, with dead-end Acyl-ACP
model = removeRxns(model, {'R01403'});

% Similar case, with dead-end Hexadecenoic acid
model = removeRxns(model, {'R08162', 'R08161'});

% Monolignol biosynthesis is cytosolic
model = removeRxns(model, {'R02381'});

% Cyanide metabolism is mitochondrial
model = removeRxns(model, {'R01931'});

% These probably only occur in bacterial metabolism
model = removeRxns(model, {'R10712', 'R11312'});

% contains 4 dead-end metabolites
model = removeRxns(model, {'R02879'});

% MEMOTE: these did not have a formula before. Taken from MetaNetX
% MNXM726595 and MNXM726611, which map to the respective MetaCyc species
model.metFormulas(strcmp(model.mets, 'C00996[h]')) = {'C34H30FeN4O4'};
model.metCharges(strcmp(model.mets, 'C00996[h]')) = -1;

model.metFormulas(strcmp(model.mets, 'C00999[h]')) = {'C34H30FeN4O4'};
model.metCharges(strcmp(model.mets, 'C00999[h]')) = -2;

% cytosolic sucrose was not actually connected to the network so far. Add
% it again with all the annotation and add a transporter
model = addCytosolMet(model, 'C00089[h]');
model.S(strcmp(model.mets, 'C00089[c]'), strcmp(model.rxns, 'Im_suc')) = 1;
model.S(strcmp(model.mets, 'Suc[c]'), strcmp(model.rxns, 'Im_suc')) = 0;
model = removeMetabolites(model, {'Suc[c]'});
model.mets(strcmp(model.mets, 'C00089[c]')) = {'Suc[c]'};
model.mets(strcmp(model.mets, 'C00089[h]')) = {'Suc[h]'};

model = addReaction(model, 'Tr_suc', 'reactionName', ...
    'Sucrose transporter', ...
    'reactionFormula', ...
    'Suc[c] --> Suc[h]', 'geneRule', 'At5g59250');
model.rxnNotes(strcmp(model.rxns, 'Tr_suc')) = ...
    {'doi.org/10.1104/pp.18.01036'};
model = addToSubsystem(model, {'Tr_suc'}, 'transport');

% The exchange is missing for this metabolite, but it has a transporter
model = addReaction(model, 'Im_oxgl', 'reactionName', ...
    'Exchange of 2-Oxoglutaramate', 'reactionFormula', '<=> C00940[c]');
model = addToSubsystem(model, {'Im_oxgl'}, 'import');
model.ub(strcmp(model.rxns, 'Im_oxgl')) = 0;

% D-Alanine is orphan so add exchange, but set bounds to 0. No evidence for
% transport however
model = addReaction(model, 'Exch_DAla', 'reactionName', ...
    'Exchange of D-alanine', 'reactionFormula', '<=> C00133[h]');
model = addToSubsystem(model, {'Exch_DAla'}, 'exchange');
model.ub(strcmp(model.rxns, 'Exch_DAla')) = 0;

% disconnected cliques; unlikely to be correct
model = removeRxns(model, {'R05409', 'R05397'});
model = removeRxns(model, {'R0488', 'R04880'});

% Not localized to chloroplast according to Uniprot Q9MB58
model = removeRxns(model, {'R13202'});

% Not localized to chloroplast according to Uniprot Q94AQ6 Q9FE17
model = removeRxns(model, {'R12391'});

% NMNH is a drug rather than a natural metabolite, according to Wikipedia
model = removeRxns(model, {'R11104'});

% Dead end even in KEGG
model = removeRxns(model, {'R04944'});

% 2-Oxosuccinamate (C02362) is produced in the peroxisome, I had added a
% transporter, but no import reaction
model = addReaction(model, 'Im_oxsuc', 'reactionName', ...
    'Import of 2-Oxosuccinamate (produced in the peroxisome)', ...
    'reactionFormula', '--> C02362[c]');
model = addToSubsystem(model, {'Im_oxsuc'}, 'import');
model.ub(strcmp(model.rxns, 'Im_oxsuc')) = 0;

%% Selenium-containing metabolites
disp('Working on selenium metabolites')
% Add import for selenate C05697 (most common form of Se) and diffusion
model = addCytosolMet(model, 'C05697[h]');
model = addReaction(model, 'Exch_SeO4', 'reactionName', ...
    'Exchange of selenate', 'reactionFormula', '<=> C05697[c]');
model = addToSubsystem(model, {'Exch_SeO4'}, 'exchange');

model = addReaction(model, 'Tr_SeO4', 'reactionName', ...
    'Diffusion of selenate (analogously to SO4)', 'reactionFormula', ...
    'C05697[c] <=> C05697[h]');
model = addToSubsystem(model, {'Tr_SeO4'}, 'transport');

% conversion to selenomethioine (R09365) is cytosolic. Therefore, remove it
% and add transport and export for its precursor selenohomocysteine C05698
model = removeRxns(model, {'R09365'});
model = addCytosolMet(model, 'C05698[h]');
model = addReaction(model, 'Tr_SeHcys', 'reactionName', ...
    'Transport of Selenohomocysteine', 'reactionFormula', ...
    'C05698[h] --> C05698[c]');
model = addToSubsystem(model, {'Tr_SeHcys'}, 'transport');
model.rxnNotes(strcmp(model.rxns, 'Tr_SeHcys')) = ...
    {'doi.org/10.1016/j.jhazmat.2020.124178'};

model = addReaction(model, 'Exch_SeHcys', 'reactionName', ...
    'Exchange of Selenohomocysteine', 'reactionFormula', '<=> C05698[c]');
model = addToSubsystem(model, {'Exch_SeHcys'}, 'exchange');
model.ub(strcmp(model.rxns, 'Exch_SeHcys')) = 0;
model.rxnNotes(strcmp(model.rxns, 'Exch_SeHcys')) = ...
    {'doi.org/10.1016/j.jhazmat.2020.124178'};

% Could not find evidence for this reaction; the enzyme included in KEGG
% seems to actually be responsible for a different reaction
model = removeRxns(model, {'R09372'});

% methylation of SeCys (C05688) to MetSeCys (C05689) is possible in
% selenium acumulating plants like broccoli
[model, addedRxnId, errorMsg] = addKEGGReactionNew(model, 'R04931', ...
    'Q4VNK0', 'doi.org/10.1016/j.jhazmat.2020.124178');
if ~isempty(errorMsg)
    rxnErrorI = rxnErrorI + 1;
    rxnsToLookAt{rxnErrorI} = [addedRxnId ': ' errorMsg];
end
model = addSubSystem(model, {'R04931'}, ...
    'Selenocompound metabolism');
model.rxnKEGGID(strcmp(model.rxns, 'R04931')) = {'R04931'};

% Gene is unknown, but it is known that conversion of Methyl-selenocysteine
% to volatile dimethyl-selenide happens in the chloroplast 
% doi.org/10.1016/j.jhazmat.2020.124178
[model, addedRxnId, errorMsg] = addKEGGReactionNew(model, 'R09369', ...
    '', 'doi.org/10.1016/j.jhazmat.2020.124178');
if ~isempty(errorMsg)
    rxnErrorI = rxnErrorI + 1;
    rxnsToLookAt{rxnErrorI} = [addedRxnId ': ' errorMsg];
end
model = addSubSystem(model, {'R09369'}, ...
    'Selenocompound metabolism');
model.rxnKEGGID(strcmp(model.rxns, 'R09369')) = {'R09369'};

% dimethyl-selenide is volatile and diffuses into the environment
model = addReaction(model, 'Ex_DMSe', 'reactionName', ...
    'Export of volatile DMSE', 'reactionFormula', 'C02535[h] -->');
model = addToSubsystem(model, {'Ex_DMSe'}, 'export');

%% Glucosamine metabolism
% R02058 is ER-localized (uniprot Q9LFU9)
% R00416 is cytosol-localized (uniprot Q940S3 O64765)
% This makes R08193 a dead-end; although it is chloroplast localized, that
% is probably because of the other reaction of the enzyme, converting 
% Glc-1-P and Glc-6-P (https://doi.org/10.3389/fpls.2024.1349064). 
% Therefore, remove it 
model = removeRxns(model, {'R02058', 'R00416', 'R08193'});

%% Choline is currently a dead-end, but it should be imported from the
% cytosol to produce glycine betaine (aka betaine): 10.1006/MBEN.2000.0158
disp('Working on choline')
choline = 'C00114[h]';
model.metCharges(strcmp(model.mets, 'C00114[h]')) = 1; % correct the charge

% Add import reaction for cytosolic choline
model = addCytosolMet(model, choline);
model = addReaction(model, 'Im_choline', 'reactionName', ...
    'Choline import', 'reactionFormula', '<=> C00114[c]');
model = addToSubsystem(model, {'Im_choline'}, 'import');

% Add choline transporter
model = addReaction(model, 'Tr_choline', 'reactionName', ...
    'Choline transporter', ...
    'reactionFormula', ...
    'C00114[c] --> C00114[h]');
model.rxnNotes(strcmp(model.rxns, 'Tr_choline')) = {'10.1006/MBEN.2000.0158'};
model = addToSubsystem(model, {'Tr_choline'}, 'transport');

% Add sink for glycine betaine since this metabolite accumulates
model = addReaction(model, 'Sk_betaine', 'reactionName', ...
    'Betaine sink to simulate accumulation', ...
    'reactionFormula', ...
    'C00719[h] --> ');
model = addToSubsystem(model, {'Sk_betaine'}, 'export');

% choline synthesis happens in the cytosol
model = removeRxns(model, {'R01030'});

% Dead ends even in KEGG
model = removeRxns(model, {'R05739', 'R00187'});

% No evidence for occurrence in the chloroplast found; also dead end in the
% maize model where this rxns came from (cpd00446[d0], cpd00697[d0])
model = removeRxns(model, {'R00191', 'R00434'});

% Incorrectly localized to chloroplast from poplar model
model = removeRxns(model, {'R00251'});

% There is no evidence for 2,5-Dioxopentanoate occurrence in chloroplasts;
% Probably a misannotation from a general reaction involving "an aldehyde"
model = removeRxns(model, {'R00264'});

% Although annotated as (among others) chloroplastic in uniprot, this is
% based on sequence prediction and not actually confirmed. Since this is a
% dead end in the chloroplast, remove it here
model = removeRxns(model, {'R00293'});

% According to KEGG, this is bidirectional
model.lb(strcmp(model.rxns, 'R00336')) = -1000;
% add a conversion between pppGpp and ppGpp, otherwise both are dead-ends,
% but there is ample evidence that they occur
[model, addedRxnId, errorMsg] = addKEGGReactionNew(model, 'R03409', '', ...
    ['GPR unknown but connects pppGpp and ppGpp, and for both there is '...
    'evidence in the chloroplast']);
if ~isempty(errorMsg)
    rxnErrorI = rxnErrorI + 1;
    rxnsToLookAt{rxnErrorI} = [addedRxnId ': ' errorMsg];
end
model.subSystems(strcmp(model.rxns, 'R03409')) = {''};
model.rxnKEGGID(strcmp(model.rxns, 'R03409')) = {'R03409'};
model.lb(strcmp(model.rxns, 'R03409')) = -1000;

% Judging by evidence, it is unlikely that R00491 occurs in the
% chloroplast, and that R00652 occurs in plants at all
model = removeRxns(model, {'R00491', 'R00652'});

% R01085 from maize model; dead end; no reviewed gene in uniprot, some
% unreviewed, most say mitochondrial
model = removeRxns(model, {'R01085'});

% Missing step in galactose metabolism; fills an important gap
[model, addedRxnId, errorMsg] = addKEGGReactionNew(model, 'R10619', ...
    '', 'Unclear enzyme, but required so connect metabolism');
if ~isempty(errorMsg)
    rxnErrorI = rxnErrorI + 1;
    rxnsToLookAt{rxnErrorI} = [addedRxnId ': ' errorMsg];
end
model = addSubSystem(model, {'R10619'}, 'Galactose metabolism');
model.rxnKEGGID(strcmp(model.rxns, 'R10619')) = {'R10619'};


%% Thiamine phosphate is chloroplastic but currently not functional
disp('Working on thiamine biosynthesis')
% thiazole part
[model, addedRxnId, errorMsg] = addKEGGReactionNew(model, 'R10711', ...
    'At5g54770', 'doi.org/10.1042/BSR20180048');
if ~isempty(errorMsg)
    rxnErrorI = rxnErrorI + 1;
    rxnsToLookAt{rxnErrorI} = [addedRxnId ': ' errorMsg];
end
model = addSubSystem(model, {'R10711'}, ...
    'MetaCyc:thiazole component of thiamine diphosphate biosynthesis III');
model = addToSubsystem(model, {'R10711'}, ...
    'Thiamine metabolism');
model.rxnKEGGID(strcmp(model.rxns, 'R10711')) = {'R10711'};

rxnIx = strcmp(model.rxns, 'R10685');
model.grRules(rxnIx) = {'At5g54770'};
model.rxnNotes(rxnIx) = {'doi.org/10.1042/BSR20180048'};
model = addToSubsystem(model, {'R10685'}, ...
    'MetaCyc:thiazole component of thiamine diphosphate biosynthesis III');
model = addToSubsystem(model, {'R10685'}, ...
    'Thiamine metabolism');
model.rxnKEGGID(rxnIx) = {'R10685'};

% TODO make sure this actually works!
% now we need a way to recover the protein-L-cysteine
% There is no such reaction in metacyc or KEGG; add the minimal balanced
% reaction
model = addReaction(model, 'protCysRec_h', 'reactionName', ...
    'unknown reaction; added to make thiamine phostphate biosynthesis possible', ...
    'reactionFormula', 'H2S[h] + C00080[h] + C21861[h] --> C02743[h]');
model = addToSubsystem(model, {'protCysRec_h'}, ...
    'MetaCyc:thiazole component of thiamine diphosphate biosynthesis III');
model = addToSubsystem(model, {'protCysRec_h'}, ...
    'Thiamine metabolism');

% pyrimidine part
model = addToSubsystem(model, {'R03472'}, ...
    'MetaCyc:4-amino-2-methyl-5-diphosphomethylpyrimidine biosynthesis I');
model.rxnKEGGID(strcmp(model.rxns, 'R03472')) = {'R03472'};
model.rxnEC(strcmp(model.rxns, 'R03472')) = {'4.1.99.17'};
model.eccodes(strcmp(model.rxns, 'R03472')) = {'4.1.99.17'};

% the reaction above produces CO and 5'-Deoxyadenosine, both of which
% cannot be metabolized at the moment
% CO has been found to be used by bacteria for reduction equivalents, but
% not in plants. Let it be exported
model = addCytosolMet(model, 'C00237[h]');
model = addReaction(model, 'Tr_CO', 'reactionName', ...
    'CO diffusion of the chloroplast', ...
    'reactionFormula', ...
    'C00237[h] --> C00237[c]');
model = addToSubsystem(model, {'Tr_CO'}, 'transport');
model.rxnNotes(strcmp(model.rxns, 'Tr_CO')) = {'Unclear fate; let CO be exported'};

model = addReaction(model, 'Ex_CO', 'reactionName', ...
    'Export of CO', 'reactionFormula', 'C00237[c] --> ');
model = addToSubsystem(model, {'Ex_CO'}, 'export');
model.lb(strcmp(model.rxns, 'Ex_CO')) = 0;

% The 5'dAdo salvage pathway likely uses the same enzymes as the methionine
% salvage pathway. This includes the following enzyme
rxnIx = strcmp(model.rxns, 'R07392');
model.grRules(rxnIx) = {[model.grRules{rxnIx} ' or At5g53850']};

% The 5'dAdo salvage pathway is not clear, but this seems to be the
% most likely route: doi.org/10.1038/s41467-018-05589-4. Not all reactions are in KEGG, however.
% First add the missing metabolites
newMets = {'C16637', 'C22280'};
for i = 1: length(newMets)
    met = newMets{i};
    res = readAndParseMet(met);
    if isempty(fieldnames(res))
        disp(['Error: no fields found for met ' met])
        errorMsg = [errorMsg 'Failed fetching metabolite details for ' met '; '];
        rxnErrorI = rxnErrorI + 1;
        rxnsToLookAt{rxnErrorI} = [addedRxnId ': ' errorMsg];
        continue
    end
    if contains(res.formula, 'n')
        errorMsg = [errorMsg met ' contains n in formula; '];
    end
    model = addMetabolite(model, [met '[h]'], res.name, ...
        res.formula, res.chebi, met, res.pubchem);
    model.metLIPIDMAPSID{end} = res.lipidmaps;
    model.metCharges(end) = 0;
end

% first rxn: phosphorylase
currRxnNote = ['Unclear reaction; our best idea so far is discussed ' ...
    'here: doi.org/10.1038/s41467-018-05589-4'];
model = addReaction(model, 'phosph_5dAdo_h', 'reactionName', ...
    '5''-Deoxyadenosine phosphorylase', ...
    'reactionFormula', 'C05198[h] + C00009[h] --> C16637[h] + C00147[h]');
model.grRules(strcmp(model.rxns, 'phosph_5dAdo_h')) = {'At5g53850'};
model.rxnNotes(strcmp(model.rxns, 'phosph_5dAdo_h')) = {currRxnNote};

% second rxn: isomerase
model = addReaction(model, 'isom_dR1P_h', 'reactionName', ...
    'D-5-Deoxyribose 1-phosphate isomerase', ...
    'reactionFormula', 'C16637[h] --> C22280[h]');
model.grRules(strcmp(model.rxns, 'isom_dR1P_h')) = {'At5g53850'};
model.rxnNotes(strcmp(model.rxns, 'isom_dR1P_h')) = {currRxnNote};

% third rxn: aldolase
model = addKEGGReactionNew(model, 'R12612', 'At5g53850', currRxnNote);

model = addSubSystem(model, {'phosph_5dAdo_h', 'isom_dR1P_h', 'R12612'}, ...
    '5''-Deoxyadenosine salvage');

% This was already there but is also part of thiamine metabolism
model = addToSubsystem(model, {'FPAL_h'}, 'Thiamine metabolism');

% the merging reaction
[model, addedRxnId, errorMsg] = addKEGGReactionNew(model, 'R10712', ...
    'At1g22940');
if ~isempty(errorMsg)
    rxnErrorI = rxnErrorI + 1;
    rxnsToLookAt{rxnErrorI} = [addedRxnId ': ' errorMsg];
end
model = addSubSystem(model, {'R10712'}, ...
    'MetaCyc:thiamine diphosphate biosynthesis IV (eukaryotes)');
model = addToSubsystem(model, {'R10712'}, ...
    'Thiamine metabolism');
model.rxnKEGGID(strcmp(model.rxns, 'R10712')) = {'R10712'};
model.eccodes(strcmp(model.rxns, 'R10712')) = {'2.5.1.3'};

% add thiamine phosphate to the precursors
model.S(strcmp(model.mets, 'C01081[h]'), ...
    strcmp(model.rxns, 'precursorPool')) = -1;

%% balancing
disp('Working on balancing')
% this balances all its reactions at once
model.metFormulas(strcmp(model.mets, 'C00154[h]')) = {'C37H62N7O17P3S'};

% tRNA ligases
% change formula of empty tRNAs to reflect the OH group that the amino acid
% binds to
model.metFormulas(endsWith(model.mets, '[h]') & ...
                  contains(model.metNames, 'tRNA(') & ...
                  ~contains(model.metNames, '-')) = repmat({'ROH'}, 20, 1);
model.metFormulas(strcmp(model.mets, 'tRNA-Glu[c]')) = {'ROH'};

% Now, the loaded tRNA formula is the amino acid plus R
model.metFormulas{strcmp(model.mets, 'C02987[h]')} = ...
	[model.metFormulas{strcmp(model.mets, 'C00025[h]')} 'R1'];
model.metFormulas{strcmp(model.mets, 'C00886[h]')} = ...
	[model.metFormulas{strcmp(model.mets, 'C00041[h]')} 'R1'];
model.metFormulas{strcmp(model.mets, 'C02163[h]')} = ...
	[model.metFormulas{strcmp(model.mets, 'C00062[h]')} 'R1'];
model.metFormulas{strcmp(model.mets, 'C03125[h]')} = ...
	[model.metFormulas{strcmp(model.mets, 'C00097[h]')} 'R1'];
model.metFormulas{strcmp(model.mets, 'C02412[h]')} = ...
	[model.metFormulas{strcmp(model.mets, 'C00037[h]')} 'R1'];
model.metFormulas{strcmp(model.mets, 'C02553[h]')} = ...
	[model.metFormulas{strcmp(model.mets, 'Ser[h]')} 'R1'];
model.metFormulas{strcmp(model.mets, 'C01931[h]')} = ...
	[model.metFormulas{strcmp(model.mets, 'Lys[h]')} 'R1'];
model.metFormulas{strcmp(model.mets, 'C02430[h]')} = ...
	[model.metFormulas{strcmp(model.mets, 'C00073[h]')} 'R1'];
model.metFormulas{strcmp(model.mets, 'C03402[h]')} = ...
	[model.metFormulas{strcmp(model.mets, 'C00152[h]')} 'R1'];
model.metFormulas{strcmp(model.mets, 'C03127[h]')} = ...
	[model.metFormulas{strcmp(model.mets, 'C00407[h]')} 'R1'];
model.metFormulas{strcmp(model.mets, 'C02554[h]')} = ...
	[model.metFormulas{strcmp(model.mets, 'C00183[h]')} 'R1'];
model.metFormulas{strcmp(model.mets, 'C02984[h]')} = ...
	[model.metFormulas{strcmp(model.mets, 'C00049[h]')} 'R1'];
model.metFormulas{strcmp(model.mets, 'C02839[h]')} = ...
	[model.metFormulas{strcmp(model.mets, 'C00082[h]')} 'R1'];
model.metFormulas{strcmp(model.mets, 'C03512[h]')} = ...
	[model.metFormulas{strcmp(model.mets, 'C00078[h]')} 'R1'];
model.metFormulas{strcmp(model.mets, 'C02702[h]')} = ...
	[model.metFormulas{strcmp(model.mets, 'C00148[h]')} 'R1'];
model.metFormulas{strcmp(model.mets, 'C06112[h]')} = ...
	[model.metFormulas{strcmp(model.mets, 'C00025[h]')} 'R1'];
model.metFormulas{strcmp(model.mets, 'C02282[h]')} = ...
	[model.metFormulas{strcmp(model.mets, 'C00064[h]')} 'R1'];
model.metFormulas{strcmp(model.mets, 'C02988[h]')} = ...
	[model.metFormulas{strcmp(model.mets, 'C00135[h]')} 'R1'];
model.metFormulas{strcmp(model.mets, 'C02047[h]')} = ...
	[model.metFormulas{strcmp(model.mets, 'C00123[h]')} 'R1'];
model.metFormulas{strcmp(model.mets, 'C03511[h]')} = ...
	[model.metFormulas{strcmp(model.mets, 'Phe[h]')} 'R1'];
model.metFormulas{strcmp(model.mets, 'C02992[h]')} = ...
	[model.metFormulas{strcmp(model.mets, 'C00188[h]')} 'R1'];

% Also balance the charges
model.metCharges(strcmp(model.mets, 'C02984[h]')) = -1;
model.metCharges(strcmp(model.mets, 'C01931[h]')) = 1;

% Use O2 and H2O to stoichiometrically balance this reaction and not
% produce NADPH in the meantime
rxnIx = strcmp(model.rxns, 'R12254');
H2OIx = find(strcmp(model.mets, 'C00001[h]'));
O2Ix = find(strcmp(model.mets, 'C00007[h]'));
model.S(H2OIx, rxnIx) = model.S(H2OIx, rxnIx) + 1;
model.S(:, rxnIx) = model.S(:, rxnIx) * 2;
model.S(O2Ix, rxnIx) = model.S(O2Ix, rxnIx) - 1;
model = addRxnNote(model, rxnIx, ['Modified to account for stoichiometry; ' ...
        'Mechanism unknown']);

% Here we have an oxygen unaccounted for, simply add molecular oxygen
% because the mechanism is unknown. Also multiply the whole reaction by 2
% to arrive at integer coefficients
rxnIx = strcmp(model.rxns, 'R08990');
model.S(:, rxnIx) = model.S(:, rxnIx) * 2;
model.S(O2Ix, rxnIx) = model.S(O2Ix, rxnIx) - 1;
model = addRxnNote(model, rxnIx, ['Modified to account for stoichiometry; ' ...
        'Mechanism unknown']);
% similar case for R13563: -2H is the balance, and the reaction proceeds
% in a direction so that these would be produced. Do not use NADHP but
% simply O2 and H2O
rxnIx = strcmp(model.rxns, 'R13563');
model.S(H2OIx, rxnIx) = model.S(H2OIx, rxnIx) + 1;
model.S(:, rxnIx) = model.S(:, rxnIx) * 2;
model.S(O2Ix, rxnIx) = model.S(O2Ix, rxnIx) - 1;
model = addRxnNote(model, rxnIx, ['Modified to account for stoichiometry; ' ...
        'Mechanism unknown']);

% similar case for R13560: -2H, +2O is the balance, and the reaction proceeds
% in a direction so that these would be produced. Do not use NADHP but
% simply O2 and H2O
rxnIx = strcmp(model.rxns, 'R13560');
model.S(H2OIx, rxnIx) = model.S(H2OIx, rxnIx) + 1;
model.S(:, rxnIx) = model.S(:, rxnIx) * 2;
model.S(O2Ix, rxnIx) = model.S(O2Ix, rxnIx) - 3;
model = addRxnNote(model, rxnIx, ['Modified to account for stoichiometry; ' ...
        'Mechanism unknown']);

% This is a dead-end even in KEGG
% For some reason this does not remove the KEGG ID entry
%rxnIx = strcmp(model.rxns, 'R03231');
%model.rxnKEGGID(rxnIx) = [];
model = removeRxns(model, {'R03231'});

% Difficult case: R08549. Follow the definition and formula of Rhea
rxnIx = strcmp(model.rxns, 'R08549');
HIx = find(strcmp(model.mets, 'C00080[h]'));
model.S(HIx, rxnIx) = 0;
model.metFormulas(strcmp(model.mets, 'C00091[h]')) = {'C25H35N7O19P3S'};

% Flavodoxin is a dead end so there is no need for it to transfer charges
% to NADP (also imbalanced reaction)
model = removeRxns(model, {'R11485'});

% First of all change bounds to produce ADP not ATP
model.ub(strcmp(model.rxns, 'R03905')) = 0;
model.lb(strcmp(model.rxns, 'R03905')) = -1000;
% Then adapt charge according to Rhea
model.metCharges(strcmp(model.mets, 'C06112[h]')) = -1;

% While these reactions do have evidence for the chloroplast in rice, no
% connection to the rest of the model does. Therefore, remove
model = removeRxns(model, {'R08363', 'R08364'});

% Unbalanced dead-end whose function can be fulfilled by R07859
model = removeRxns(model, {'R07860'});

% Similar to other 15-cis-phytoene desaturase (1.3.5.5) reactions, use
% plastoquinone to deliver the H
rxnIx = strcmp(model.rxns, 'R07510');
plastoquinolIx = find(strcmp(model.mets, 'C16695[h]'));
plastoquinoneIx = find(strcmp(model.mets, 'C10385[h]'));
model.S(plastoquinolIx, rxnIx) = model.S(plastoquinolIx, rxnIx) + 2;
model.S(plastoquinoneIx, rxnIx) = model.S(plastoquinoneIx, rxnIx) - 2;
model = addRxnNote(model, rxnIx, ['Modified to account for stoichiometry; ' ...
        'Mechanism unknown']);

% Adapt the formulas to reflect Rhea and balance the reactions
model.metFormulas(strcmp(model.mets, 'C00681[h]')) = {'C4H6O7PR'};
model.metFormulas(strcmp(model.mets, 'C00416[h]')) = {'C5H5O8PR2'};
model.metFormulas(strcmp(model.mets, 'C04899[h]')) = {'C21H28N7O18P3'};
model.metFormulas(strcmp(model.mets, 'C04752[h]')) = {'C6H8N3O7P2'};
model.metFormulas(strcmp(model.mets, 'C00054[h]')) = {'C10H11N5O10P2'};

model.metCharges(strcmp(model.mets, 'NH4[c]')) = 1;

% This one has too many H already, so the balancing is fine
rxnIx = strcmp(model.rxns, 'R01777');
model.S(HIx, rxnIx) = model.S(HIx, rxnIx) - 5;

rxnIx = strcmp(model.rxns, 'R12172');
model.S(HIx, rxnIx) = model.S(HIx, rxnIx) - 4;

% Balance reactions with protons and water
model = automaticBalancing(model);

%% miscellaneous
disp('Working on miscellaneous changes')
% Let MEMOTE find the NGAM reaction
model.rxns(strcmp(model.rxns, 'R00086')) = {'NGAM_h'};

% check for error messages
if ~isempty(rxnsToLookAt)
    disp(rxnsToLookAt)
end

%% Unblocking of reactions
disp('Unblocking reactions')

% naringenin is not produced in the chloroplast 
model = removeRxns(model, {'R02446', 'R01613'});

% connects alpha-zeacarotene with delta-carotene
[model, addedRxnId, errorMsg] = addKEGGReactionNew(model, 'R06961', '', ...
    'GPR unknown but connects the carotenoid metabolism');
if ~isempty(errorMsg)
    rxnErrorI = rxnErrorI + 1;
    rxnsToLookAt{rxnErrorI} = [addedRxnId ': ' errorMsg];
end
model = addToSubsystem(model, {'R06961'}, 'Carotenoid metabolism');
model.rxnKEGGID(strcmp(model.rxns, 'R06961')) = {'R06961'};
% But: is unbalanced by 2 H; balance using water and oxygen to prevent
% production of reducing equivalents
rxnIx = strcmp(model.rxns, 'R06961');
model.S(H2OIx, rxnIx) = model.S(H2OIx, rxnIx) + 1;
model.S(:, rxnIx) = model.S(:, rxnIx) * 2;
model.S(O2Ix, rxnIx) = model.S(O2Ix, rxnIx) - 1;

% remove reactions that contain more than one dead-end metabolites since it
% is unlikely that they can be connected
while true
    bool_S = (model.S ~= 0);
    disconnectedMetsIx = (sum(bool_S, 2) < 2);
    numDisconnectedMets = sum(bool_S(disconnectedMetsIx, :));
    if sum(numDisconnectedMets > 1)
        model = removeRxns(model, model.rxns(numDisconnectedMets > 1));
    else
        break
    end
end

% TODO check all addKEGG and addReaction if I update the rxn2subs 
% TODO make pathways work with the SBML model
% TODO add all new import and exchange reactions to the media tsv
%% write pcm to disk
disp('Writing to disk')
if isfield(model, 'A')
    model = rmfield(model, 'A');
end
if isfield(model, 'C')
    model = rmfield(model, 'C');
end
model = updateFromGrRules(model);

writeCbModel(model, 'fileName', ...
    [modelPath 'pcm.v' num2str(importVersion + 1) '.mat']);
writeCbModel(model, 'fileName', ...
    [modelPath 'pcm.v' num2str(importVersion + 1) '.xml']);

%% Functions
function comp = extractComp(metId)
    metComp = regexp(metId, '\[(.*?)\]', 'tokens');
    if ~isempty(metComp)
        comp = metComp{end}{1};
    else
        comp = '';
    end
end

function model = addRxnNote(model, rxnIx, rxnNote)
    if ~isempty(model.rxnNotes{rxnIx})
        model.rxnNotes{rxnIx} = [model.rxnNotes{rxnIx} '; ' rxnNote];
    else
        model.rxnNotes{rxnIx} = rxnNote;
    end
end


function res = readAndParseMet(metID)

    res = struct();

    % webrequest
    [response, success] = ...
        tryRequest('https://rest.kegg.jp/get/compound:', metID);
    if ~success
        return
    end

    % parsing
    lines = strsplit(response, '\n');

    res.name = '';
    res.formula = '';
    res.pubchem = '';
    res.chebi = '';
    res.lipidmaps = '';

    i = 0;
    while i < length(lines)
        i = i + 1;
        line = lines{i};
        
        % Extract the name
        if startsWith(line, 'NAME')
            nameLine = line(length('NAME') + 1:end);

            % if name is not terminated, check next lines
            while ~contains(nameLine, ';')
                if i+1 > length(lines) || ~startsWith(lines{i+1}, ' ')
                    break
                end
                i = i + 1;
                nameLine = [nameLine, ' ', strtrim(lines{i})];
            end

            % use the string collected so far
            if contains(nameLine, ';')
                res.name = strtrim(nameLine(1:strfind(nameLine, ';')-1));
            else
                res.name = strtrim(nameLine);
            end
        end
    
        % Extract the formula
        if startsWith(line, 'FORMULA')
            res.formula = strtrim(line(length('FORMULA') + 1:end));
        end
    
        % Extract the DBLINKS values
        if startsWith(line, 'DBLINKS')
            dbLinksLine = strtrim(line(length('DBLINKS') + 1:end));

            % collect whole dblinks string
            while true
                if i+1 > length(lines) || ~startsWith(lines{i+1}, ' ')
                    break
                else
                    i = i + 1;
                    dbLinksLine = [dbLinksLine, '$$$', strtrim(lines{i})];
                end
            end

            % check parts
            dbLinksParts = strsplit(dbLinksLine, '$$$');
            for j = 1:length(dbLinksParts)
                part = strtrim(dbLinksParts{j});
                if startsWith(part, 'PubChem:')
                    res.pubchem = strtrim(part(strfind(part, ':')+1:end));
                elseif startsWith(part, 'ChEBI:')
                    tmp = strtrim(part(strfind(part, ':')+1:end));
                    chebis = strsplit(tmp, ' ');
                    res.chebi = chebis{1};
                elseif startsWith(part, 'LipidMaps:')
                    res.lipidmaps = strtrim(part(strfind(part, ':')+1:end));
                end
            end
        end
    end
end


function [res, success] = tryRequest(base_url, ele)
    tries = 0;
    res = '';
    while true
        try
            tries = tries + 1;
            res = webread([base_url ele]);
            success = true;
            return
        catch ME
            if strcmp(ME.identifier, 'MATLAB:webservices:HTTP403StatusCodeError') && tries < 5
                pause(15);
            else
                disp(['Error during webread of ' ele ' after ' ...
                      num2str(tries) ' tries'])
                success = false;
                return
            end
        end
    end
end