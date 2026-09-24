clearvars -except gurobiAvailable projDir; clc;

modelPath = strcat(projDir, 'Data/pcm/');
importVersion = 1;
load([modelPath 'pcm.v' num2str(importVersion) '.mat'], 'model');

%% General changes
disp('Working on general changes')
% Update the subsystem of import reactions
model = addToSubsystem(model, {'Im_suc', 'Im_str', 'Im_o2'}, 'import');

% Make the upper bound of these import non-zero to ensure flux is possible
% through every reaction in default medium
model.ub(strcmp(model.rxns, 'Im_suc')) = 0.1;
model.ub(strcmp(model.rxns, 'Im_str')) = 0.1;
model.ub(strcmp(model.rxns, 'Im_o2')) = 0.1;

% G6P (beta-G6P) annotation with C00092 (general G6P) leads to problems in
% plugAndPlay because C00092 exists as its own metabolite
model.metKEGGID{strcmp(model.mets, 'G6P[h]')} = 'C01172';

% Add source for a bunch of transporters that I had curated:
model.rxnNotes(strcmp(model.rxns, 'NDH_h')) = {'10.7554/eLife.49305'};
model.rxnNotes(strcmp(model.rxns, 'Tr_xanth')) = {'10.1073/pnas.2502160122'};
model.rxnNotes(strcmp(model.rxns, 'Tr_THF')) = {'10.1111/j.1399-3054.2006.00587.x'};
model.rxnNotes(strcmp(model.rxns, 'Tr_pABA')) = {'10.1111/j.1399-3054.2006.00587.x'};

% CTP import had a low upper bound to prevent the cell from using CTP for
% energy metabolism. From this paper however 
% (doi.org/10.1105/tpc.112.096743) it is apparent that in
% order for pyrimidines to be produced in the cytosol and mitochondria,
% carbamoyl-aspartate has to be exported from the plastid. While the
% cytidine transport route has not been identified, this gives reasoning to
% couple CAs export to Ura import. PRPP is also needed but the carbon atoms
% are not retained in uracil, so do not couple it. Since the route of
% cytidine import is unknown, keep the solution of low upper bounds
model = addMetabolite(model, 'CONST_pyri[h]', ...
    'Pyrimidine (precursor) ex-/import pseudo metabolite');
model.S(strcmp(model.mets, 'CONST_pyri[h]'), ...
    strcmp(model.rxns, 'Tr_Cas_h')) = 1;
model.S(strcmp(model.mets, 'CONST_pyri[h]'), ...
    strcmp(model.rxns, 'Tr_Ura')) = -1;
model.rxnNotes(strcmp(model.rxns, 'Tr_Ura')) = {'doi.org/10.1105/tpc.112.096743'};
% pyrimidine steps not in the chloroplast
model = removeRxns(model, {'R01993', 'R01870', 'R00965'});

% After extensive research on pyrimidine metabolism in chloroplasts, I
% found that there are single enzymes that are known, but their connections
% are missing. Generally, it is believed that cytidine metabolites are
% esupplied from the cytosolic biosynthesis. Therefore, add transporters
% for the dead-end cytidine metabolites instead of assuming chloroplast
% reactions
model = addImportTransport(model, 'C00881', 'Deoxycytidine', 'dCyt', ...
    ['Cytidine synthesis occurs in the cytosol exclusively, but the ' ...
    'reactions downstream are identified in chloroplast']);
model.ub(strcmp(model.rxns, 'Im_dCyt')) = 0.1;
model = addImportTransport(model, 'C00239', 'dCMP', 'dCMP', ...
    ['Cytidine synthesis occurs in the cytosol exclusively, but the ' ...
    'reactions downstream are identified in chloroplast']);
model.ub(strcmp(model.rxns, 'Im_dCMP')) = 0.1;

% These import reactions are set by the medium, but also set their upper
% bound to less in the general model
model.ub(strcmp(model.rxns, 'Im_GDP')) = 0.1;
model.ub(strcmp(model.rxns, 'Im_GTP')) = 0.1;
model.ub(strcmp(model.rxns, 'Exch_C00008[c]')) = 0.1;
model.ub(strcmp(model.rxns, 'Im_AMP')) = 0.1;

%% Changes to address dead-ends
disp('Addressing dead-ends')
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
                           'R00857', 'R08618', 'R01340', 'R01341', ...
                           'R07256'});

% remove metabolites that are outside any 0 reactions (these are in the
% cytosol)
model = removeMetabolites(model, {'C20694[c]', 'C12287[c]'});

% reaction from xenobiotics and carcinogenesis
model = removeRxns(model, {'R07113'});

% To be able to use FAD, let it interchange reduction equivalents with NAD
% and NADP
[model, addedRxnId, errorMsg] = addKEGGReactionNew(model, 'R09520', 'unknown', ...
    'Added for modeling purpose. Some way of activating FAD must be present');
if ~isempty(errorMsg)
    rxnErrorI = rxnErrorI + 1;
    rxnsToLookAt{rxnErrorI} = [addedRxnId ': ' errorMsg];
end
model = addToSubsystem(model, {'R09520'}, 'Riboflavin metabolism');
model.lb(strcmp(model.rxns, 'R09520')) = -1000;
model.rxnKEGGID(strcmp(model.rxns, 'R09520')) = {'R09520'};

[model, addedRxnId, errorMsg] = addKEGGReactionNew(model, 'R09748', 'unknown', ...
    'Added for modeling purpose. Some way of activating FAD must be present');
if ~isempty(errorMsg)
    rxnErrorI = rxnErrorI + 1;
    rxnsToLookAt{rxnErrorI} = [addedRxnId ': ' errorMsg];
end
model = addToSubsystem(model, {'R09520'}, 'Riboflavin metabolism');
model.lb(strcmp(model.rxns, 'R09748')) = -1000;
model.rxnKEGGID(strcmp(model.rxns, 'R09748')) = {'R09748'};

% The siroheme biosynthesis pathway's occurrence in plants is supported by 
% metacyc and occurrence of the following reactions in the chloroplast is 
% supported by Uniprot. But it is a protein, so it is out of scope for this
% model
% #Siroheme
model = removeRxns(model, {'R02864', 'R03194'});

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
model.ub(strcmp(model.rxns, 'Im_oxgl')) = 0.1;

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
model.ub(strcmp(model.rxns, 'Im_oxsuc')) = 0.1;

% beta-alanine biosynthesis is relatively well characterized
% (doi.org/10.3389/fpls.2019.00921), but its subcellular localization I had
% to gather from uniprot entries
% Enzymes producing 3-Aminopropanal are annotated in peroximsome and
% cytosol, but the enzyme converting it to beta-alanine is plastidic.
% Therefore, add a 3-Aminopropanal import and transporter
model = addCytosolMet(model, 'C05665[h]');
model = addReaction(model, 'Im_3APr', 'reactionName', ...
    '3-Aminopropanal import', 'reactionFormula', '<=> C05665[c]');
model = addToSubsystem(model, {'Im_3APr'}, 'import');

model = addReaction(model, 'Tr_3APr', 'reactionName', ...
    '3-Aminopropanal transporter', ...
    'reactionFormula', ...
    'C05665[c] --> C05665[h]');
model.rxnNotes(strcmp(model.rxns, 'Tr_3APr')) = {['Uniprot ' ...
    'information of subcellular localization of enzymes producing ' ...
    'or consuming 3-Aminopropanal']};
model = addToSubsystem(model, {'Tr_3APr'}, 'transport');

% Found no convincing evidence
model = removeRxns(model, {'R00919', 'R00995', 'R01150', 'R01224'});

% Evidence (uniprot) points to peroxisome
model = removeRxns(model, {'R00927'});

% Since evidence exists for a reaction producing D-glutamate (At5g57850),
% add the racemase as the solution of least-assumptions. There are some
% plant proteins annotated with that enzymatic function
[model, addedRxnId, errorMsg] = addKEGGReactionNew(model, 'R00260', ...
    'A0A9Q0BYR0', ['Evidence exists only for production of ' ...
    'D-glutamate, so there must be some way to get rid of it']);
if ~isempty(errorMsg)
    rxnErrorI = rxnErrorI + 1;
    rxnsToLookAt{rxnErrorI} = [addedRxnId ': ' errorMsg];
end
model.rxnKEGGID(strcmp(model.rxns, 'R00260')) = {'R00260'};

% The annotated genes were different methyltransferases. The correct
% protein At2g32160 has no evidence for chloroplast-localization
model = removeRxns(model, {'R01159'});

% No evidence for chloroplast localization
model = removeRxns(model, {'R01167', 'R01718', 'R01678', 'R02147'});

% At1g31190 is the correct gene, impl2 is another reaction. But, this is
% likely a side-reaction of the enzyme that is not relevant since the
% substrate is probably not present in chloroplasts
model = removeRxns(model, {'R01185'});

% From maize model, but uniprots knows no protein in maize, and in other
% plants this is localized to the nucleus (since this is a signaling
% molecule)
model = removeRxns(model, {'R01232'});

% Evidence in Maize points to cytoplasm and secretion
model = removeRxns(model, {'R01378', 'R01372'});

% evidence points to ER
model = removeRxns(model, {'R01281'});

% no evidence for chloroplast found
model = removeRxns(model, {'R01210', 'R01355', 'R01736'});

% No evidence for mannose metabolism/import into chloroplasts and this is a
% dead end
model = removeRxns(model, {'R01329', 'R01326'});

% No evidence found; proteins rather point to other compartments
model = removeRxns(model, {'R01364', 'R01360'});

% Evidence rather points towards cytosol and peroxisome
model = removeRxns(model, {'R01422', 'R02453'});

% These are connected with cell wall biosynthesis, which is not connected
% to the chloroplast
model = removeRxns(model, {'R01384', 'R01473'});

% Found no evidence for this reactions (also pyrimidine)
model = removeRxns(model, {'R01567'});

% These are vacuolar / ER / cytoplasmic / peroxisomal
model = removeRxns(model, {'R02976', 'R02464', 'R01468', 'R01274'});

% Monoloignol biosynthesis is cytosolic, specialized pathways have been
% found in peroxisome and chloroplasts, but not these basic ones
model = removeRxns(model, {'R01616', 'R01617', 'R03337', 'R03339'});

% Rather cytosolic (uniprot Q949W8)
model = removeRxns(model, {'R01639'});

% There is no evidence for plants
model = removeRxns(model, {'R01645'});

% These are the remaining sequential reactions of the identified enzyme
% K7WQ45 that produce nerylneryl diphosphate 
% (doi.org/10.1371/journal.pone.0119302). This is a precursor for
% many non-canonical diterpenes (doi.org/10.1186/s12870-020-2293-x). 
% Because these are species-specific and not too well studied, kepp the
% pathway around but end at that precursor
[model, addedRxnId, errorMsg] = addKEGGReactionNew(model, 'R12638', ...
    'K7WQ45', 'doi.org/10.1371/journal.pone.0119302');
if ~isempty(errorMsg)
    rxnErrorI = rxnErrorI + 1;
    rxnsToLookAt{rxnErrorI} = [addedRxnId ': ' errorMsg];
end
model.rxnKEGGID(strcmp(model.rxns, 'R12638')) = {'R12638'};

[model, addedRxnId, errorMsg] = addKEGGReactionNew(model, 'R12639', ...
    'K7WQ45', 'doi.org/10.1371/journal.pone.0119302');
if ~isempty(errorMsg)
    rxnErrorI = rxnErrorI + 1;
    rxnsToLookAt{rxnErrorI} = [addedRxnId ': ' errorMsg];
end
model.subSystems(strcmp(model.rxns, 'R12639')) = {''};
model.rxnKEGGID(strcmp(model.rxns, 'R12639')) = {'R12639'};

model = addReaction(model, 'Ex_NNPP', 'reactionName', ...
    'Export of nerylneryl-diphosphate to keep the pathway', ...
    'reactionFormula', 'C21784[h] -->');
model = addToSubsystem(model, {'Ex_NNPP'}, 'export');
model.rxnNotes(strcmp(model.rxns, 'Ex_NNPP')) = {['Non-canonical precursor of' ...
    ' highly species-specific diterpenes ' ...
    '(doi.org/10.1186/s12870-020-2293-x)']};

% Difficult case: Q8W593 gives evidence for de-toxification of
% methylglyoxal to S-lactoylglutathione (SLG). The canonical route of
% disposal of SLG would be cleavage into glutathione (GSH) and D-lactate.
% Alternatively, SLG might be exported into other compartments, like
% mitochondria where that cleavage has been confirmed. However, there is no
% evidence for either option. GSH is needed in chloroplasts for redox
% balance, so direct cleavage would make sense. There is no evidence for a
% further metabolism of D-lactate in chloroplasts, which would make the
% export of SLG more likely. Because including chloroplastic cleavage would
% require even more assumptions, I decide for the simple export of SLG.
model = addExportTransport(model, 'C03451', 'Lactoylglutathione', 'SLG', ...
    ['Evidence for production exists (Q8W593); further pathway unclear; ' ...
    'export for metabolism in mitochondria']);
% remove the cleavage reaction
model = removeRxns(model, {'R01736'});
% add an import for methylglyoxal so de-toxification can run
model = addImportTransport(model, 'C00546', 'Methylglyoxal', 'MetGlyox', ...
    'There is evidence for de-toxification of Methylglyoxal, so add import');

% No evidence for this reaction
model = removeRxns(model, {'R01817'});
% This would leave some reactions around mannose-6P blocked; R00772 is
% actually bidirectional according to MetaCyc MANNPISOM-RXN and unblocks
model.lb(strcmp(model.rxns, 'R00772')) = -1000;

% remove At3g20540; is DNA-polymerase?
% remove At1g27680; transfers adenyl instead of uridinyl
% At3g56040 is the only sensible gene; transfers UDP to glucose and R02634
% might be side reaction
model.grRules(strcmp(model.rxns, 'R02634')) = {'At3g56040'};
% However, the connecting reactions result in dead end and have no real
% support; were included from poplar model and no literature evidence was
% found. Therefore, remove the whole pathway
model = removeRxns(model, {'R01980', 'R02634', 'R01385', 'R00286'});

% R02013, limonene production, is legit but is not metabolized but stored
% or evaporates. Add export
model = addReaction(model, 'Ex_lim', 'reactionName', ...
    'Export of limonene', 'reactionFormula', 'C00521[h] -->');
model = addToSubsystem(model, {'Ex_lim'}, 'export');

% Not in chloroplast
model = removeRxns(model, {'R02082'});

% Dead end and no evidence
model = removeRxns(model, {'R02089'});

% trans-cinnamate is produced in the cytosol
model = removeRxns(model, {'R00697', 'R02256', 'R02254'});

% Found no evidence this clique around phenylacetaldehyde occurs in
% chloroplasts
model = removeRxns(model, {'R02536', 'R02611', 'R01377', 'R02613'});

% Found no sufficient evidence for these reactions
model = removeRxns(model, {'R03347', 'R07411', 'R07419', 'R02741', ...
                           'R02425', 'R07273', 'R05706', 'R08658', ...
                           'R10049', 'R09121', 'R03171', 'R02529', ...
                           'R05615', 'R03546', 'R05616', 'R02323', ...
                           'R02719', 'R08086', 'R07405', 'R10054', ...
                           'R10563', 'R06954'});

% Trehalose (also called mycose) is indeed produced in different plants, 
% and used for various functions. Add a transport and export.
model = addExportTransport(model, 'C01083', 'alpha,alpha-Trehalose', ...
            'mycose', '');

% This reaction is legit andcan be connected with the same enzyme
% A0AAV9CR38 that is already present in the model (but not annotated)
pcm.grRules(strcmp(model.rxns, 'R03067')) = {'A0AAV9CR38'};

% This part of folate biosynthesis is cytosolic. The annotated protein
% A0A699GEZ9 does not seem to actually support the reactions
model = removeRxns(model, {'R03503', 'R03067'});
% This makes R02237 a dead-end reaction. Because it is only one of the
% reactions catalyzed by At5g05980, remove it
model = removeRxns(model, {'R02237'});

% Lienar pathway as part of branched-chain amino acid degradation, which is
% unlikely to be correctly placed in chloroplasts. The best evidence is
% provided by A0A9E7H647, an unreviewed uniprot entry. Remove the whole
% pathway
model = removeRxns(model, {'R03869', 'R05066', 'R05064', 'R04224'});
% This makes R03869 a dead-end reaction. Because it is only one of the
% reactions catalyzed by At4g34240, remove it as well
model = removeRxns(model, {'R03869'});

% GDP-fucose de novo synthesis. There is some hints at a chloroplast
% localization (Q9SNY3, A0A7G2DUG9). Add the respective genes and an export
% function for GDP fucose, since it is required for cell wall biosynthesis
% (doi.org/10.1104/pp.103.022368)
model.grRules(strcmp(model.rxns, 'R00888')) = ...
    {[model.grRules{strcmp(model.rxns, 'R00888')} ' or At5g66280']};
model = addExportTransport(model, 'C00325', 'GDP-L-fucose', 'GDPfuc', ...
    'required for cell wall biosynthesis doi.org/10.1104/pp.103.022368');

% These might not be per se wrong, but are disconnected and the only
% curated gene associated to these () also support other reactions.
% According to MetaCyc, the reaction removed here is not considered to be
% physiologically relevant
% (https://biocyc.org/reaction?orgid=META&id=R145-RXN)
model = removeRxns(model, {'R07392', 'R07393', 'R07394'});

% The two reactions R03815 and R07618 seem to be the same
metAnnoFields = {'metSmiles', 'metInChIString', 'metKEGGID', ...
                 'metChEBIID', 'metPubChemID', 'metMetaNetXID', ...
                 'metBioCycID', 'metLIPIDMAPSID'};
model = mergeMets(model, 'C15973[h]', 'DHL[h]', metAnnoFields);
model = mergeMets(model, 'C15972[h]', 'C02051[h]', metAnnoFields);
model = removeRxns(model, {'R07618'});

% Likely correct, but I see no possible route for its product octanoic acid
% (C06423). It might be involved in coniine biosynthesis, but subcellular
% localization of this compound is unclear doi.org/10.1111/febs.13410
% Therefore, simply export
% Same for R08158, decanoic acid
% Same for R04014, dodecanoic acid
% Same for R08159, tetradecanoic acid
model = addExportTransport(model, 'C06423', 'Octanoic acid', 'octa', ...
    'Evidence for enzymes producing is there, but furhter fate is unclear');
model = addExportTransport(model, 'C01571', 'Decanoic acid', 'deca', ...
    'Evidence for enzymes producing is there, but furhter fate is unclear');
model = addExportTransport(model, 'C02679', 'Dodecanoic acid', 'dodeca', ...
    'Evidence for enzymes producing is there, but furhter fate is unclear');
model = addExportTransport(model, 'C06424', 'Tetradecanoic acid', 'tetradeca', ...
    'Evidence for enzymes producing is there, but furhter fate is unclear');

% 9beta-Pimara-7,15-diene is a precursor for momilactones, but its fate is
% not well known, and so is the subcellular localization of subsequent 
% enzmes. Add an export function
model = addExportTransport(model, 'C18225', '9beta-Pimara-7,15-diene', ...
    'pimDiene', ...
    ['precursor for oryzallexins, but further pathway''s sub-cellular ' ...
    'localization unknown doi.org/10.1111/j.1365-313x.2010.04408.x']);
% Same for Sandaracopimaradiene and oryzallexins
model = addExportTransport(model, 'C11877', 'Sandaracopimaradiene', ...
    'sandar', ...
    ['precursor for oryzallexins, but further pathway''s sub-cellular ' ...
    'localization unknown dx.doi.org/10.1104/pp.111.187518']);
% Same for ent-Isokaurene
model = addExportTransport(model, 'C20145', 'ent-Isokaurene', ...
    'entIkaur', ...
    ['Diterpene synthesized in chloroplasts, but further pathway''s' ...
    ' sub-cellular localization unknown 10.1016/j.febslet.2011.09.038']);

%  remove because dAMP is a dead-end with no evidence for reactions (take
%  place in mitochondria instead
model = removeRxns(model, {'R01547'});

% dead-end side reactions of the associated enzyme(s)
model = removeRxns(model, {'R02527', 'R02528', 'R03096', 'R02678' ...
        'R03921', 'R05052'});

% associated enzmes canonically work on Fru1,6BP instead of Fru1P. Dead-end
% side-reaction
model = removeRxns(model, {'R02568'});

% Merge beta- with the general fructose
model = mergeMets(model, 'C02336[h]', 'C00095[h]', metAnnoFields);

% No evidence for the chloroplast and dead-end
model = removeRxns(model, {'R03293', 'R03291', 'R00471'});

% dead-end; reaction definition of the enzyme is extremely broad (same 
% enzyme(s) for all of these reactions); likely involved in de-toxification
% of herbicides and not relevant for the normal metabolism, which is the 
% intended use of the pcm
model = removeRxns(model, {'R07003', 'R07004', 'R07023', 'R07024', ...
    'R07025', 'R07026', 'R07092', 'R07094', 'R07116'});

% enzyme not yet characterized
[model, addedRxnId, errorMsg] = addKEGGReactionNew(model, 'R07508', ...
    'unknown', ['Enzyme not yet characterized but beta-tocopherol can ' ...
    'be produced and alpha-tocopherol is the most active vitamin E']);
if ~isempty(errorMsg)
    rxnErrorI = rxnErrorI + 1;
    rxnsToLookAt{rxnErrorI} = [addedRxnId ': ' errorMsg];
end
model = addToSubsystem(model, {'R07508'}, ...
    'MetaCyc:vitamin E biosynthesis (tocopherols)');
model.rxnKEGGID(strcmp(model.rxns, 'R07508')) = {'R07508'};
% The reaction is unbalanced; Assume SAM (aMet) as methyldonor
rxnIx = strcmp(model.rxns, 'R07508');
model.S(strcmp(model.mets, 'aMet[h]'), rxnIx) = ...
    model.S(strcmp(model.mets, 'aMet[h]'), rxnIx) - 1;
model.S(strcmp(model.mets, 'C00080[h]'), rxnIx) = ...
    model.S(strcmp(model.mets, 'C00080[h]'), rxnIx) + 1;
model.S(strcmp(model.mets, 'aH-Cys[h]'), rxnIx) = ...
    model.S(strcmp(model.mets, 'aH-Cys[h]'), rxnIx) + 1;

% Found some more general reactions. Concrete instances are present except
% for TDP, which does not exist in the model and would be dead end
model = removeRxns(model, {'R07261_2', 'R00331', 'R02320', 'R07261', ...
    'R00331'});

% 2-trans-rest-cis-prenyl-diphosphates are only seen in bacteria. Plant
% equivalent is all-trans-prenyl-diphosphates, so this is likely an error
% from homology assignment. 
model = removeRxns(model, {'R07269', 'R08753', 'R08752', 'R08751', ...
    'R08750', 'R08749', 'R08748', 'R05555'});

% This made me find some wrong gprs:
model.grRules(strcmp(model.rxns, 'FPPS_h')) = ...
    {strrep(model.grRules{strcmp(model.rxns, 'FPPS_h')}, 'At4g38460', ...
    '( At4g38460 and At4g36810)')};
model.grRules(strcmp(model.rxns, 'GPPS_h')) = ...
    {strrep(model.grRules{strcmp(model.rxns, 'GPPS_h')}, 'At4g38460', ...
    '( At4g38460 and At4g36810)')};

% from R07267 R05613 R05611 R05612, remove At5g19040 At1g68460 At3g63110
% because these catalyze different reactions
genes = {'At5g19040', 'At1g68460', 'At3g63110'};
rxns = {'R07267', 'R05613', 'R05611', 'R05612'};
for rIx = 1:length(rxns)
    currRxnIx = strcmp(model.rxns, rxns{rIx});
    for gIx = 1:length(genes)
        model.grRules{currRxnIx} = ...
            strrep(model.grRules{currRxnIx}, [genes{gIx} ' or '], '');
        model.grRules{currRxnIx} = ...
            strrep(model.grRules{currRxnIx}, [' or ' genes{gIx}], '');
    end
end

% Astaxanthin is used for pigmentation of flowers etc, not relevant for
% light harvesting in leafs https://doi.org/10.1093/plphys/kiab428. It is
% produced in green algae, but not in chloroplasts 
% https://doi.org/10.1111/tpj.12713
model = removeRxns(model, {'R07568', 'R07572'});

% The reaction is a dead end in KEGG, and uniprot rather hints at a
% function of the associated gene in modifying proteins
model = removeRxns(model, {'R07606'});

% Linalool is a volatile terpene and aromatic compound. It is produced by
% many plants in the chloroplasts. Its derivatives in wikipedia are not in
% KEGG and the derivatives in KEGG are not supported for plants. Add an
% export
model = addReaction(model, 'Ex_Linalool', 'reactionName', ...
    'Export of volatile linalool, also simulating accumulation', ...
    'reactionFormula', 'C11389[h] -->');
model = addToSubsystem(model, {'Ex_Linalool'}, 'export');

% R07674 is a dead-end side reactions and R07678 has no support for
% chloroplasts
model = removeRxns(model, {'R07674', 'R07678'});

% None of the associated genes are correct. Betalains occur in some plants
% of the order Caryophyllales but I could not find genes and no chloroplast
% localization
model = removeRxns(model, {'R08836'});

% two-sided dead end cliques without support for chloroplast
model = removeRxns(model, {'R02265', 'R02268'});
model = removeRxns(model, {'R02258', 'R02261'});

% two-sided dead end clique; CoA biosynthesis occurs in the cytosol
model = removeRxns(model, {'R02472', 'R02473', 'R03018'});
% CoA biosynthesis occurs in the cytosol; No evidence for chloroplast
% localization of these reactions; dead-end clique
model = removeRxns(model, {'R03269', 'R03035', 'R00130'});

% two-sided dead end clique; related to cell wall biosynthesis; chloroplast
% localization unlikely
model = removeRxns(model, {'R06514', 'R02777'});

% Found no evidence for chloroplast occurrence
model = removeRxns(model, {'R02926', 'R01787'});

% Fate of beta-alanine in chloroplasts is unclear, but it is known to
% accumulate for stress resistance. Add export to reflect that
model = addReaction(model, 'Ex_bAla', 'reactionName', ...
    'Export of beta alanine to reflect accumulation', ...
    'reactionFormula', 'C00099[h] -->');
model.rxnNotes(strcmp(model.rxns, 'Ex_bAla')) = ...
    {['known to accumulate for stress resistance' ...
    'https://doi.org/10.3389/fpls.2019.00921']};
model = addToSubsystem(model, {'Ex_bAla'}, 'export');

% slim evidence for chloroplast localization and no way was found to
% produce inositol-1,4-P
model = removeRxns(model, {'R03393', 'R01186'});

% Both are part of Yang cycle for methioine recycling
% R01401 no evidence for chloroplast; rather cytsol. R04143 no known 
% subcellular localization, but occurrence was  found in vascular tissue
% and thus not relevant for out leaf-chloroplast model. R00179 (ACS) was
% found predominantly in the cytosol.
% https://doi.org/10.48130/ph-0025-0007 
model = removeRxns(model, {'R01401', 'R04143', 'R00179'});

% This makes R00997 dead-end, but it is most likely mitochondrial anyway
% dx.doi.org/10.1111/j.1742-4658.2005.04567.x
model = removeRxns(model, {'R00997'});

% Degradation of cell wall constituent (not modelled) and therefore
% dead-end
model = removeRxns(model, {'R01101', 'R05549'});

%% Degradation of metabolites produced spontaneously or outside the model
disp(['Working on degradation of various metabolites ' ...
    'produced spontaneously or outside the model'])

% D-Alanine can be degraded, but is an orphan metabolite. Since these 
% metabolite might occur randomly, add exchange but set bounds to 0 to
% prevent the model from using it for energy metabolism. No evidence for 
% transport found
model = addReaction(model, 'Exch_DAla', 'reactionName', ...
    'Exchange of D-alanine', 'reactionFormula', '<=> C00133[h]');
model = addToSubsystem(model, {'Exch_DAla'}, 'exchange');
model.ub(strcmp(model.rxns, 'Exch_DAla')) = 0;
% No evidence found for degradation reactions of D-arginine and analogous
% reaction with L-alanine
model = removeRxns(model, {'R02923', 'R08197'});

% pseudouridine is a similar case: It is formed post-transcriptionally in
% RNA, but during degradation, it has to be recycled (which it is).
% Modification of RNA is outside the scope of the model, so add an exchange
% reaction
model = addReaction(model, 'Exch_pseudouridine', 'reactionName', ...
    ['Exchange of Pseudouridine from degradation of post-' ...
    'transcriptionally modified RNA'], 'reactionFormula', '<=> C02067[h]');
model = addToSubsystem(model, {'Exch_pseudouridine'}, 'exchange');
model.ub(strcmp(model.rxns, 'Exch_pseudouridine')) = 0.1;

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

% Detoxification reactions of hydrogen cyanide
model = addReaction(model, 'Exch_HCN', 'reactionName', ...
    'Exchange of HCN', 'reactionFormula', '<=> C01326[h]');
model = addToSubsystem(model, {'Exch_HCN'}, 'exchange');
model.ub(strcmp(model.rxns, 'Exch_HCN')) = 0.1;
model = addSubSystem(model, {'R03524', 'R01267', 'Exch_HCN'}, ...
    'HCN detoxification');

%% Chlorophyll degradation
% It progresses until "primary fluorescent chlorophyll catabolite" in 
% plastids and then is exported to the ER/vacuole
% doi.org/10.1093/pcp/pcae093
disp('Working on chlorophyll degradation')

[model, addedRxnId, errorMsg] = addKEGGReactionNew(model, 'R09033', ...
    'At4g11910', 'doi.org/10.1093/pcp/pcae093');
if ~isempty(errorMsg)
    rxnErrorI = rxnErrorI + 1;
    rxnsToLookAt{rxnErrorI} = [addedRxnId ': ' errorMsg];
end
model.rxnKEGGID(strcmp(model.rxns, 'R09033')) = {'R09033'};

[model, addedRxnId, errorMsg] = addKEGGReactionNew(model, 'R09032', ...
    'At4g37000', 'doi.org/10.1093/pcp/pcae093');
if ~isempty(errorMsg)
    rxnErrorI = rxnErrorI + 1;
    rxnsToLookAt{rxnErrorI} = [addedRxnId ': ' errorMsg];
end
model.rxnKEGGID(strcmp(model.rxns, 'R09032')) = {'R09032'};

% add the spontaneous reaction R11228 to fix the dead-end R11227 / C18022
[model, addedRxnId, errorMsg] = addKEGGReactionNew(model, 'R11228', ...
    'spontaneous');
if ~isempty(errorMsg)
    rxnErrorI = rxnErrorI + 1;
    rxnsToLookAt{rxnErrorI} = [addedRxnId ': ' errorMsg];
end
model.rxnKEGGID(strcmp(model.rxns, 'R11228')) = {'R11228'};

model = addExportTransport(model, 'C18098', ['primary fluorescent' ...
    ' chlorophyll catabolite'], 'chlCat', ...
    'doi.org/10.1093/pcp/pcae093');

% add a subsystem for chlorophyll degradation
model = addSubSystem(model, {'R09032', 'R09033', 'R08921', 'R11228', ...
    'Tr_chlCat'}, 'chlorophyll degradation');

%% Updates to the biomass
disp('Updating biomass')

% dTTP is a dead end currently; add nucleotides to the biomass
model = addMetabolite(model, 'pNc[h]', 'Nucleotides lump metabolite');
model = addReaction(model, 'nucleotidePool', 'reactionName', ...
    'nucelotidePool (biomass)', 'reactionFormula', ...
    ['0.31 C00131[h] + ' ... % dATP
    '0.19 C00458[h] +' ... % dCTP
    ' 0.18 C00286[h] + ' ... % dGTP
    '0.16 C00459[h] + ' ... % dTTP
    '0.16 C00460[h] ' ... % dUTP
    '<=> pNc[h]']);
model = addToSubsystem(model, {'nucleotidePool'}, 'pseudo reaction');

model.S(strcmp(model.mets, 'pNc[h]'), ...
        strcmp(model.rxns, 'BiomassRxn')) = -1;

% add another pigment to the pigment biomass: (3Z)-Phytochromobilin
model.S(strcmp(model.mets, 'C05913[h]'), ...
        strcmp(model.rxns, 'PigmentPool')) = -0.1;

% epsilon carotene is a pigment and has been found in multiple plants,
% including zea mays
% https://pubchem.ncbi.nlm.nih.gov/compound/epsilon-Carotene
model.S(strcmp(model.mets, 'C16276[h]'), ...
        strcmp(model.rxns, 'PigmentPool')) = -0.1;

% 7,8-dihydro-beta-carotene is the endpoint of a pigment biosynthesis
% pathway that is present in many species, including A. thaliana:
% https://doi.org/10.1105/tpc.8.9.1613
model.S(strcmp(model.mets, 'C16291[h]'), ...
        strcmp(model.rxns, 'PigmentPool')) = -0.1;

% S-sulfo-L-cysteine is an important signaling molecule
% https://doi.org/10.1093/mp/sst168
model.S(strcmp(model.mets, 'C05824[h]'), ...
        strcmp(model.rxns, 'precursorPool')) = -0.1;

% ent-kaurene is the GA precursor, and the downstream steps are catalyzed
% by an enzme facing outside the chloroplast 
% doi.org/10.1046/j.1365-313X.2001.01150.x
model = addExportTransport(model, 'C06090', 'ent-kaurene', 'kaur', ...
    'doi.org/10.1046/j.1365-313X.2001.01150.x');
model.S(strcmp(model.mets, 'C06090[c]'), ...
        strcmp(model.rxns, 'precursorPool')) = -0.1;

% vitamin E is synthesized and predominantly located in the chloroplast.
% Although the transporter outside of the chloroplast has not been
% identified, it has been found that transport must occur:
% doi.org/10.1016/j.tplants.2019.08.006
model = addExportTransport(model, 'C02477', 'alpha-Tocopherol', 'vitE', ...
    'doi.org/10.1016/j.tplants.2019.08.006');
model.S(strcmp(model.mets, 'C02477[c]'), ...
        strcmp(model.rxns, 'precursorPool')) = -0.1;

% Also tocotrienols (other form of vitamin E) are produced in chloroplasts
% doi.org/10.1146/annurev.arplant.56.032604.144301
rids_5_5_1_24 = {'R10624', 'R10623'};
for rIx = 1:length(rids_5_5_1_24)
    rid = rids_5_5_1_24{rIx};
    [model, addedRxnId, errorMsg] = addKEGGReactionNew(model, rid, ...
        'At4g32770', 'doi.org/10.1146/annurev.arplant.56.032604.144301');
    if ~isempty(errorMsg)
        rxnErrorI = rxnErrorI + 1;
        rxnsToLookAt{rxnErrorI} = [addedRxnId ': ' errorMsg];
    end
    model.rxnKEGGID(strcmp(model.rxns, rid)) = {rid};
end
model = addToSubsystem(model, rids_5_5_1_24, ...
    'MetaCyc:vitamin E biosynthesis (tocotrienols)');

rids_2_1_1_95 = {'R10491', 'R10492'};
for rIx = 1:length(rids_2_1_1_95)
    rid = rids_2_1_1_95{rIx};
    [model, addedRxnId, errorMsg] = addKEGGReactionNew(model, rid, ...
        'At1g64970', 'doi.org/10.1146/annurev.arplant.56.032604.144301');
    if ~isempty(errorMsg)
        rxnErrorI = rxnErrorI + 1;
        rxnsToLookAt{rxnErrorI} = [addedRxnId ': ' errorMsg];
    end
    model.rxnKEGGID(strcmp(model.rxns, rid)) = {rid};
end
model = addToSubsystem(model, rids_2_1_1_95, ...
    'MetaCyc:vitamin E biosynthesis (tocotrienols)');
% For this one, add an export to simulate accumulation, but do not add to
% biomass, since tocotrienol forms mainly in nonphotosynthetic tissues
% doi.org/10.1016/j.tplants.2019.08.006
model = addReaction(model, 'Ex_aTocotrienol', 'reactionName', ...
    'Export of alpha tocotrienol to simulate accumulation', ...
    'reactionFormula', 'C14153[h] -->');
model.rxnNotes(strcmp(model.rxns, 'Ex_aTocotrienol')) = ...
    {'doi.org/10.1016/j.tplants.2019.08.006'};
model = addToSubsystem(model, {'Ex_aTocotrienol'}, 'export');

model = addReaction(model, 'Ex_bTocotrienol', 'reactionName', ...
    'Export of beta tocotrienol to simulate accumulation', ...
    'reactionFormula', 'C14154[h] -->');
model.rxnNotes(strcmp(model.rxns, 'Ex_bTocotrienol')) = ...
    {'doi.org/10.1016/j.tplants.2019.08.006'};
model = addToSubsystem(model, {'Ex_bTocotrienol'}, 'export');

% PLP is synthesized in the chloroplast and helps with resistance to
% photooxidative damage https://doi.org/10.1016/j.plaphy.2010.10.003
% But, biosynthesis (R10089) happens in the cytosol
% (https://doi.org/10.1073/pnas.0506228102) and downstream processes are 
% catalyzed in other organelles as far as I found
% doi.org/10.1093/plcell/koae176 
% Therefore, import it, ensure interconversion of B6 vitamers works by 
% adding B6 to biomass as required for normal chloroplast functionality. 
model = removeRxns(model, {'R10089'});
% Gprs have unclear connection to the reactions. This pathway is bacterial
model = removeRxns(model, {'R05086', 'R05085'});

model = addImportTransport(model, 'C00018', ...
    'Pyridoxal phosphate; vitamin B6', 'PLP', ['Biosynthesis is ' ...
    'cytosolic (doi.org/10.1073/pnas.0506228102) but B6 and ' ...
    'enzymes interconverting its vitamers have been found in ' ...
    'chloroplasts (e.g. doi.org/10.1016/j.plaphy.2010.10.003)']);

model = addMetabolite(model, 'B6[h]', 'metName', 'B6 vitamers');
model.S(strcmp(model.mets, 'B6[h]'), ...
    strcmp(model.rxns, 'precursorPool')) =  -1;
model = addReaction(model, 'b6Pool', 'reactionName', ...
    'Pool for B6 vitamers', 'reactionFormula', ...
    'C00534[h] + C00647[h] + C00250[h] + C00018[h] + C00314[h] + C00627[h]  <=> B6[h]', ...
    'geneRule', 'PSEUDO');
model = addToSubsystem(model, {'b6Pool'}, 'pseudo reaction');

[model, addedRxnId, errorMsg] = addKEGGReactionNew(model, 'R00173', ...
    'At2g33255', 'doi.org/10.1093/plphys/kiac048');
if ~isempty(errorMsg)
    rxnErrorI = rxnErrorI + 1;
    rxnsToLookAt{rxnErrorI} = [addedRxnId ': ' errorMsg];
end
model.rxnKEGGID(strcmp(model.rxns, 'R00173')) = {'R00173'};
model = addToSubsystem(model, {'R00173'}, ...
    'Vitamin B6 metabolism');

[model, addedRxnId, errorMsg] = addKEGGReactionNew(model, 'R01710', ...
    'At5g49970', '');
if ~isempty(errorMsg)
    rxnErrorI = rxnErrorI + 1;
    rxnsToLookAt{rxnErrorI} = [addedRxnId ': ' errorMsg];
end
model.rxnKEGGID(strcmp(model.rxns, 'R01710')) = {'R01710'};
model = addToSubsystem(model, {'R01710'}, ...
    'Vitamin B6 metabolism');
model.lb(strcmp(model.rxns, 'R01710')) = -1000;

%% Selenium-containing metabolites
disp('Working on selenium metabolites')
% Add import for selenate C05697 (most common form of Se) and diffusion
model = addCytosolMet(model, 'C05697[h]');
model = addReaction(model, 'Exch_SeO4', 'reactionName', ...
    'Exchange of selenate', 'reactionFormula', '<=> C05697[c]');
model = addToSubsystem(model, {'Exch_SeO4'}, 'exchange');
model = addToSubsystem(model, {'Exch_SeO4'}, 'Selenocompound metabolism');

model = addReaction(model, 'Tr_SeO4', 'reactionName', ...
    'Diffusion of selenate (analogously to SO4)', 'reactionFormula', ...
    'C05697[c] <=> C05697[h]');
model = addToSubsystem(model, {'Tr_SeO4'}, 'transport');
model = addToSubsystem(model, {'Tr_SeO4'}, 'Selenocompound metabolism');

% conversion to selenomethioine (R09365) is cytosolic. Therefore, remove it
% and add transport and export for its precursor selenohomocysteine C05698
model = removeRxns(model, {'R09365'});
model = addCytosolMet(model, 'C05698[h]');
model = addReaction(model, 'Tr_SeHcys', 'reactionName', ...
    'Transport of Selenohomocysteine', 'reactionFormula', ...
    'C05698[h] --> C05698[c]');
model = addToSubsystem(model, {'Tr_SeHcys'}, 'transport');
model = addToSubsystem(model, {'Tr_SeHcys'}, 'Selenocompound metabolism');
model.rxnNotes(strcmp(model.rxns, 'Tr_SeHcys')) = ...
    {'doi.org/10.1016/j.jhazmat.2020.124178'};

model = addReaction(model, 'Ex_SeHcys', 'reactionName', ...
    'Export of Selenohomocysteine', 'reactionFormula', 'C05698[c] -->');
model = addToSubsystem(model, {'Ex_SeHcys'}, 'exchange');
model = addToSubsystem(model, {'Ex_SeHcys'}, 'Selenocompound metabolism');
model.ub(strcmp(model.rxns, 'Ex_SeHcys')) = 0.1;
model.rxnNotes(strcmp(model.rxns, 'Ex_SeHcys')) = ...
    {'doi.org/10.1016/j.jhazmat.2020.124178'};

% No evidence for production of cystathionine (or selenocystathionine) from
% acetyl-homoserine in chloroplasts (or plants) found.
model = removeRxns(model, {'R04945', 'R03217'});

% Now, C05699 Selenocystathionine cannot be produced, so add reaction that
% uses Phosphohomoserine for synthesis, analogously to cystathionine
% synthesis
model = addReaction(model, 'SeCTHS_h', 'reactionName', ...
    'plant cystathionine gamma-synthase', ...
    'reactionFormula','C05688[h] + C01102[h] --> C00009[h] + C05699[h]',...
    'geneRule', 'At3g01120');
model.rxnEC{strcmp(model.rxns, 'SeCTHS_h')} = '2.5.1.48';
model.eccodes(strcmp(model.rxns, 'R10712')) = {'2.5.1.48'};
model = addToSubsystem(model, {'SeCTHS_h'}, 'Selenocompound metabolism');

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
model.rxnKEGGID(strcmp(model.rxns, 'R09369')) = {'R09369'};
model = addToSubsystem(model, {'R09369', 'R04931'}, ...
    'Selenocompound metabolism');

% dimethyl-selenide is volatile and diffuses into the environment
model = addReaction(model, 'Ex_DMSe', 'reactionName', ...
    'Export of volatile DMSE', 'reactionFormula', 'C02535[h] -->');
model = addToSubsystem(model, {'Ex_DMSe'}, 'export');
model = addToSubsystem(model, {'Ex_DMSe'}, 'Selenocompound metabolism');

%% Glucosamine metabolism
% R02058 is ER-localized (uniprot Q9LFU9)
% R00416 is cytosol-localized (uniprot Q940S3 O64765)
% This makes R08193 a dead-end; although it is chloroplast localized, that
% is probably because of the other reaction of the enzyme, converting 
% Glc-1-P and Glc-6-P (https://doi.org/10.3389/fpls.2024.1349064). 
% Therefore, remove it 
% R01965 R01961 are dead-end side-reactions of hexokinases
% R00768 is now a dead-end and no evidence for chloroplast was found
model = removeRxns(model, {'R02058', 'R00416', 'R08193', 'R01965', ...
    'R01961', 'R00768'});

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
    'unknown', 'Unclear enzyme, but required to connect metabolism');
if ~isempty(errorMsg)
    rxnErrorI = rxnErrorI + 1;
    rxnsToLookAt{rxnErrorI} = [addedRxnId ': ' errorMsg];
end
model = addToSubsystem(model, {'R10619'}, 'Galactose metabolism');
model.rxnKEGGID(strcmp(model.rxns, 'R10619')) = {'R10619'};

% These are detoxification reactions (or from microbial/animal metabolism)
% that are of no relevance in natural chloroplast metabolism
model = removeRxns(model, {'R00703', 'R00754'});


%% Melatonin biosynthesis
disp('Working on melatonin biosynthesis')
% tryptophan biosynthesis occurs in chloroplast, already in the model.
% Conversion to tryptamine (R00685) occurs in the cytosol, as does the
% conversion of N-Acetylserotonin (NAS) to melatonin (R03130) and
% conversion of serotonin to 5-Methoxytryptamine (R02912) (remove later to
% keep the metabolites for now)

% Then, the step of SNAT: serotonin to NAS or 5-methoxytryptamine to
% melatonin, R02911 or [not in KEGG], occurs in the chloroplast again.
model = addToSubsystem(model, {'R02911'}, ...
    'MetaCyc:serotonin and melatonin biosynthesis II');

model = addReaction(model, 'SNAT_Mel_h', 'reactionName', ...
    'SNAT reaction producing melatonin', ...
    'reactionFormula', ...
    'C05659[h] + C00024[h] --> C01598[h] + CoA[h] + C00080[h]', ...
    'geneRule', 'At1g32070');
model = addToSubsystem(model, {'SNAT_Mel_h'}, ...
    'MetaCyc:serotonin and melatonin biosynthesis II');
model.rxnNotes(strcmp(model.rxns, 'SNAT_Mel_h')) = ...
    {'https://doi.org/10.1111/jpi.12364'};

% Therefore, link the import of serotonin or 5-methoxytryptamine
% to the export of tryptophan
model = addMetabolite(model, 'CONST_trp[h]', ...
    'Melatonin biosynthesis ex-/import constraint');
model = addCytosolMet(model, 'C00780[h]'); % serotonin
model = addCytosolMet(model, 'C05659[h]'); % 5-methoxytryptamine

model = addReaction(model, 'Tr_serotonin', 'reactionName', ...
    'transport of serotonin', ...
    'reactionFormula', ...
    'C00780[c] + CONST_trp[h] --> C00780[h]');
model = addToSubsystem(model, {'Tr_serotonin'}, 'transport');
model = addToSubsystem(model, {'Tr_serotonin'}, ...
    'MetaCyc:serotonin and melatonin biosynthesis II');
model.rxnNotes(strcmp(model.rxns, 'Tr_serotonin')) = ...
    {'https://doi.org/10.1111/jpi.12364'};

model = addReaction(model, 'Tr_5MT', 'reactionName', ...
    'transport of 5-methoxytryptamine', ...
    'reactionFormula', ...
    'C05659[c] + CONST_trp[h] --> C05659[h]');
model = addToSubsystem(model, {'Tr_5MT'}, 'transport');
model = addToSubsystem(model, {'Tr_5MT'}, ...
    'MetaCyc:serotonin and melatonin biosynthesis II');
model.rxnNotes(strcmp(model.rxns, 'Tr_5MT')) = ...
    {'https://doi.org/10.1111/jpi.12364'};

% Although this is not physiological (because the chloroplast provides all
% tryptophan of the cell), require tryptophan export and N-Acetylserotonin
% import to happen in a 1:1 ratio, since these constraints will be removed
% during plugAndPlay.m anyway
model = addCytosolMet(model, 'C00078[h]');
model = addReaction(model, 'Tr_Trp', 'reactionName', ...
    'transport of tryptophan', ...
    'reactionFormula', ...
    'C00078[h] --> C00078[c] + CONST_trp[h]');
model = addToSubsystem(model, {'Tr_Trp'}, 'transport');
model = addToSubsystem(model, {'Tr_Trp'}, ...
    'MetaCyc:serotonin and melatonin biosynthesis II');
model.rxnNotes(strcmp(model.rxns, 'Tr_Trp')) = ...
    {'https://doi.org/10.1111/jpi.12364'};

% Now add the exchange reactions
model = addReaction(model, 'Ex_Trp', 'reactionName', ...
    'Export of tryptophan', 'reactionFormula', 'C00078[c] -->');
model = addToSubsystem(model, {'Ex_Trp'}, 'export');

model = addReaction(model, 'Im_5MT', 'reactionName', ...
    'Import of 5-methoxytryptamine', 'reactionFormula', '--> C05659[c]');
model = addToSubsystem(model, {'Im_5MT'}, 'import');
model = addReaction(model, 'Im_serotonin', 'reactionName', ...
    'Import of serotonin', 'reactionFormula', '--> C00780[c]');
model = addToSubsystem(model, {'Im_serotonin'}, 'import');

% add transport and export of products out of the chloroplast
% use one pseudo metabolite to indicate that melatonin or its precursor was
% exported
model = addMetabolite(model, 'CONST_mel_prec[h]', ...
    'Melatonin precursor export pseudo metabolite');
model = addCytosolMet(model, 'C00978[h]');
model = addReaction(model, 'Tr_NAS', 'reactionName', ...
    'transport of N-Acetylserotonin', ...
    'reactionFormula', ...
    'C00978[h] --> C00978[c] + CONST_mel_prec[h]');
model = addToSubsystem(model, {'Tr_NAS'}, 'transport');
model = addToSubsystem(model, {'Tr_NAS'}, ...
    'MetaCyc:serotonin and melatonin biosynthesis II');
model.rxnNotes(strcmp(model.rxns, 'Tr_NAS')) = ...
    {'https://doi.org/10.1111/jpi.12364'};

model = addReaction(model, 'Ex_NAS', 'reactionName', ...
    'Export of N-Acetylserotonin', 'reactionFormula', 'C00978[c] -->');
model = addToSubsystem(model, {'Ex_NAS'}, 'export');

model = addCytosolMet(model, 'C01598[h]');
model = addReaction(model, 'Tr_melatonin', 'reactionName', ...
    'transport of melatonin', ...
    'reactionFormula', ...
    'C01598[h] --> C01598[c] + CONST_mel_prec[h]');
model = addToSubsystem(model, {'Tr_melatonin'}, 'transport');
model = addToSubsystem(model, {'Tr_melatonin'}, ...
    'MetaCyc:serotonin and melatonin biosynthesis II');
model.rxnNotes(strcmp(model.rxns, 'Tr_melatonin')) = ...
    {'https://doi.org/10.1111/jpi.12364'};

model = addReaction(model, 'Ex_melatonin', 'reactionName', ...
    'Export of melatonin', 'reactionFormula', 'C01598[c] -->');
model = addToSubsystem(model, {'Ex_melatonin'}, 'export');

% Finally, given melatonin importance in signalling and defense 
% (https://doi.org/10.1111/jpi.12364), add export of it or its precursor to 
% demand reaction 
model.S(strcmp(model.mets, 'CONST_mel_prec[h]'), ...
    strcmp(model.rxns, 'precursorPool')) = -0.1;

% Now remove the cytosolic reactions
model = removeRxns(model, {'R00685', 'R03130', 'R02912'});

%% Thiamine-P biosynthesis is chloroplastic but currently not functional
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
model = addToSubsystem(model, {'R10712'}, ...
    'MetaCyc:thiamine diphosphate biosynthesis IV (eukaryotes)');
model = addToSubsystem(model, {'R10712'}, ...
    'Thiamine metabolism');
model.rxnKEGGID(strcmp(model.rxns, 'R10712')) = {'R10712'};
model.rxnEC(strcmp(model.rxns, 'R10712')) = {'2.5.1.3'};
model.eccodes(strcmp(model.rxns, 'R10712')) = {'2.5.1.3'};

% add thiamine phosphate to the precursors
model.S(strcmp(model.mets, 'C01081[h]'), ...
    strcmp(model.rxns, 'precursorPool')) = -0.1;

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

%% miscellaneous
disp('Working on miscellaneous changes')
% Let MEMOTE find the NGAM reaction
model.rxns(strcmp(model.rxns, 'R00086')) = {'NGAM_h'};

% according to metacyc and kegg, this should proceed towards dNDP
model.lb(strcmp(model.rxns, 'R02017')) = -1000;
model.ub(strcmp(model.rxns, 'R02017')) = 0;
model.lb(strcmp(model.rxns, 'R02018')) = -1000;
model.ub(strcmp(model.rxns, 'R02018')) = 0;
model.lb(strcmp(model.rxns, 'R02019')) = -1000;
model.ub(strcmp(model.rxns, 'R02019')) = 0;
model.lb(strcmp(model.rxns, 'R02024')) = -1000;
model.ub(strcmp(model.rxns, 'R02024')) = 0;

% according to MetaCyc and KEGG Pathways, these should be reversible
model.lb(strcmp(model.rxns, 'R08676')) = -1000;
model.ub(strcmp(model.rxns, 'R08676')) = 1000;
model.lb(strcmp(model.rxns, 'R00806')) = -1000;
model.ub(strcmp(model.rxns, 'R00806')) = 1000;

% some misannotated genes


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

[model, addedRxnId, errorMsg] = addKEGGReactionNew(model, 'R01920', 'At1g23820', ...
    ['Subcellular localization might be chloroplastic but is not ' ...
    'confirmed. The rxn completes the pathway around putrescine, so ' ...
    'its occurrence seems likely']);
if ~isempty(errorMsg)
    rxnErrorI = rxnErrorI + 1;
    rxnsToLookAt{rxnErrorI} = [addedRxnId ': ' errorMsg];
end
model.lb(strcmp(model.rxns, 'R01920')) = -1000;
model.rxnKEGGID(strcmp(model.rxns, 'R01920')) = {'R01920'};

% Putrescine can only be produced so far
% ODC, producing putrescine from ornithine, is in the ER:
%   doi.org/10.1016/j.plantsci.2024.112232
model = removeRxns(model, {'R00670'});

% But the synthesis from arginine is localized to the chloroplast.
% Downstream enzymes (NATA1, SPDS1/2) are in other compartments however.
 model = addExportTransport(model, 'C00134', 'Putrescine', 'putr', ...
    'Necessary based on enzymes producing/consuming putrescine');

% there is evidence for an ODC in chloroplasts in some plants (  ), but
% downstream metabolism is not clear whatsoever. Also, generally, its
% localization is thought to be ER doi.org/10.1016/j.plantsci.2024.112232.
% Therefore, remove these reactions as unlikely to be present in most
% chloroplasts and poorly understood pathways
model = removeRxns(model, {'R00462', 'R06740'});

% Function of At3g47450 (R00557 is the same as R00111+R00558) as NOS has 
% been questioned doi.org/10.1074/jbc.M804838200 and NO is currently a 
% dead-end in the model.
% Although NO is known to be produced in chloroplasts and be involved in 
% regulation, signaling, photodamage resilience 
% (doi.org/10.1093/jxb/eraa504), its production mechanism is unknown.
% Therefore, remove all. If kept, should remove lump rxn R00557 and make
% the other two reversible 
model = removeRxns(model, {'R00557', 'R00111', 'R00558'});

% No evidence for even plants
model = removeRxns(model, {'R00123'});

% Long, blocked pathway with reactions from branched chain amino acid
% degradation. All rxns are not chloroplastic:
% R02085 probably mitochondrial
% R00238 cytosol, peroxisome
% R01978 Unclear loc but probably mitochondrion
% R04138 gpr probably misannotation and rxn rather mitochondrial
% R04095 mitochondrial
% R01651 No evidence even for plants
model = removeRxns(model, {'R02085', 'R00238', 'R01978', 'R04138', ...
    'R04095', 'R01651'});

% Rather mitochondrial
model = removeRxns(model, {'R01868'}); 

% Succinate can be produced but cannot be converted further. In higher
% plants, fate is unclear: Either conversion to malate, or export to
% cytosol. In Chlamydomonas, chloroplastic succinate dehydrogenase (Suc ->
% Fum) has been confirmed doi.org/10.1104/pp.90.3.1084 
% This might not hold for land plants, and no enzymes are known. Therefore,
% export succinate, and also keep rxns without gpr in this union-model
model = addExportTransport(model, 'C00042', 'Succinate', 'Succ', ...
    'Fate of succinate in chloroplasts of plants is unknown; export it');
model.grRules(strcmp(model.rxns, 'R02164')) = {'unknown'};

% This part of photorespiration is clearly peroxisomal!
model = removeRxns(model, {'R00475', 'R00717'});
% allantoate metabolism is not chloroplastic: peroxisome and ER
model = removeRxns(model, {'R02422', 'R05554', 'R02423', 'R00469'});

% Identical to R03050, but corresponding to KEGG definition
model = removeRxns(model, {'R04672'});

% Enzyme not identified in plants, but is required for thiamine metabolism   
[model, addedRxnId, errorMsg] = addKEGGReactionNew(model, 'R00615', 'unknown', ...
    'Enzyme unknown. Rxn required for thiamine metabolism');
if ~isempty(errorMsg)
    rxnErrorI = rxnErrorI + 1;
    rxnsToLookAt{rxnErrorI} = [addedRxnId ': ' errorMsg];
end
model = addToSubsystem(model, {'R00615'}, 'Thiamine metabolism');
model.rxnKEGGID(strcmp(model.rxns, 'R00615')) = {'R00615'};

% No evidence for this whole clique; integration was based on an automatic
% annotation of 2.4.1.82 to 4 compartments, one of which was chloroplast.
model = removeRxns(model, {'R02411', 'R03418', 'R01103', ...
    'R03634', 'R01194'});

% Found no evidence for these, and I do not see a way how its product,
% hypoxanthine, can not be a dead-end but connected to the metabolism
model = removeRxns(model, {'R01560', 'R01770', 'R01863', 'R01126'});

% R07511 is a lump of R09658+R09656, and all are reversible. Remove the
% lump to prevent futile cycle
model = removeRxns(model, {'R07511'});

% Multiple reactions are blocked because the pair quinone / hydroquinone
% can only be reduced, coming back to quinone is impossible. Proceed
% analogously to the PLM with plastoquinone and use O2 and At4g22260
model = addReaction(model, 'AOX4_quin_h', 'reactionName', ...
    'Alternative NAD(P)H-ubiquinone oxidoreductase', ...
    'reactionFormula', 'C00007[h] + 2 C15603[h] --> 2 C00001[h] + 2 C15602[h]');
model.grRules(strcmp(model.rxns, 'AOX4_quin_h')) = {'At4g22260'};
model.subSystems(strcmp(model.rxns, 'AOX4_quin_h')) = {''};

% pairs of dead-end, side reactions of At4g34240 (promiscuos enzyme, 
% moonlighting activity)
model = removeRxns(model, {'R02549', 'R01986'});
model = removeRxns(model, {'R04904', 'R04903'});

% pairs of dead-end, side reactions of multiple enzymes
model = removeRxns(model, {'R05395', 'R05396'});

% triangle clique (R08869+R08868 = R05165) that is completely disconnected.
% Has reviewed uniprot entry but no publication associated and localization
% is from a bulk-MS study
model = removeRxns(model, {'R08869', 'R08868', 'R05165'});

% two-sided dead-end cliques with no evidence
model = removeRxns(model, {'R02124', 'R08379'});
model = removeRxns(model, {'R03581', 'R10308'});
model = removeRxns(model, {'R03544', 'R03545'});
model = removeRxns(model, {'R07141', 'R07142'});

% dead-end, very slim evidence for chloroplast localization
model = removeRxns(model, {'R02300', 'R02301'});

% two-sided dead-end, all 4 are the same reaction
model = removeRxns(model, {'R03239', 'R03237', 'R03236', 'R03238'});

% Add a piece of evidence for myo-inositol production in chloroplasts
model.rxnNotes(strcmp(model.rxns, 'R07324')) = ...
    {[model.rxnNotes{strcmp(model.rxns, 'R07324')} ...
    'Enzyme unknown, but presence in plant chloroplasts known:' ...
    ' doi.org/10.1111/j.1365-3040.1996.tb00023.x']};
% Main production of myo-inositol is cytosolic, but there are also
% chloroplast enzymes known, and myo-inositol has been found to be shared
% across the whole plant doi.org/10.1105/tpc.10.5.753 . Export it 
model = addExportTransport(model, 'C00137', 'myo-Inositol', 'MI', ...
    ['Known to be shared in the whole plant doi.org/10.1105/tpc.10.5.753' ...
    ' doi.org/10.1111/j.1365-3040.1996.tb00023.x']);

% Heptaprenyl-diphosphate is a valid intermediate compound, but no KEGG
% reaction exists for its further processing into nonaprenyl-diphosphate 
% (or first octa, which can be processed into nona), which is required for
% plastoquinol-9 biosynthesis
model = addReaction(model, 'hepta2octaPP_h', 'reactionName', ...
    ['(E)-heptaprenyl-diphosphate:isopentenyl-diphosphate ' ...
    'heptaprenyltranstransferase'], ...
    'reactionFormula', 'C00129[h] + C04216[h] --> C00013[h] + C04146[h]');
model.subSystems(strcmp(model.rxns, 'hepta2octaPP_h')) = {''};
model.grRules(strcmp(model.rxns, 'hepta2octaPP_h')) = ...
    model.grRules(strcmp(model.rxns, 'R07267'));

%% Address stoichiometric consistency
% Actually, there is nothing to do here since the stoichiometric
% inconsistency comes only from pseudo-metabolites and photons, which is
% fine

%% Final, automatic procedures that have to be at the end
disp('Working on final, automatic procedures')

% Update metComps array
for i = 1:length(model.mets)
    model.metComps(i) = find(strcmp(model.comps, extractComp(model.mets{i})));
end

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

% Balance reactions with protons and water
model = automaticBalancing(model);

% Fill KEGG rxn IDs for all rxns that have KEGG IDs as their IDs
for rxnIx = 1:length(model.rxns)
    rid = model.rxns{rxnIx};
    
    if isempty(model.rxnKEGGID{rxnIx})
        % Determine if rxn ID is a KEGG ID
        if ~isempty(regexp(rid, '^R\d+$', 'once'))
            model.rxnKEGGID{rxnIx} = rid;
        end
    end
end

% TODO add all new import and exchange reactions to the media tsv
% TODO remove stoichiometrically identical rxns

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

% For xml: Change format of the subsystems array
% This will give warnings with writeCbModel, but will map subsystems
% correctly in the .xml
xmlModel = model;
for i = 1:length(xmlModel.rxns)
    s = xmlModel.subSystemNames(logical(xmlModel.rxn2subSystem(i, :)));
    xmlModel.subSystems{i} = s;
end

% Also name the EC field so COBRA detects it
xmlModel.rxnECNumbers = xmlModel.rxnEC;

writeCbModel(xmlModel, 'fileName', ...
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

function model = addImportTransport(model, metId, metName, shortName, ...
    trRxnNote)
    if ~exist('trRxnNote','var')
        trRxnNote=char.empty;
    end

    model = addCytosolMet(model, [metId '[h]']);
    model = addReaction(model, ['Tr_' shortName], 'reactionName', ...
        ['Transport of ' metName ' from rest of the cell'], ...
        'reactionFormula', ...
        [metId '[c] --> ' metId '[h]']);
    model = addToSubsystem(model, {['Tr_' shortName]}, 'transport');
    model.rxnNotes(strcmp(model.rxns, ['Tr_' shortName])) = {trRxnNote};
    
    model = addReaction(model, ['Im_' shortName], 'reactionName', ...
        ['Import of ' metName], 'reactionFormula', ['--> ' metId '[c]']);
    model = addToSubsystem(model, {['Im_' shortName]}, 'import');
end

function model = addExportTransport(model, metId, metName, shortName, ...
    trRxnNote)
    if ~exist('trRxnNote','var')
        trRxnNote=char.empty;
    end

    model = addCytosolMet(model, [metId '[h]']);
    model = addReaction(model, ['Tr_' shortName], 'reactionName', ...
        ['Transport of ' metName ' to rest of the cell'], ...
        'reactionFormula', ...
        [metId '[h] --> ' metId '[c]']);
    model = addToSubsystem(model, {['Tr_' shortName]}, 'transport');
    model.rxnNotes(strcmp(model.rxns, ['Tr_' shortName])) = {trRxnNote};
    
    model = addReaction(model, ['Ex_' shortName], 'reactionName', ...
        ['Export of ' metName], 'reactionFormula', [metId '[c] -->']);
    model = addToSubsystem(model, {['Ex_' shortName]}, 'export');
end