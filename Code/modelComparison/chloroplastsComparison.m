% gather data about the different models used for comparison
clearvars -except gurobiAvailable projDir; clc;

pcmImportVersion = 1;

%% import models for comparison
aragem = readCbModel([projDir 'Data/comparison_models/' 'AraGEM_valid.xml']);
maize = readCbModel([projDir 'Data/comparison_models/' 'maize.valid.xml']);
rice = readCbModel([projDir 'Data/comparison_models/' 'iOS2164.sbml']);
mintzOron = readCbModel([projDir 'Data/comparison_models/' 'mintz-oron.leaf_model.valid.sbml']);
panAlgae = readCbModel([projDir 'Data/comparison_models/' 'panAlgae.wKEGG.xml']);
poplar = readCbModel([projDir 'Data/comparison_models/' 'poplar.cobraOut.xml']);
potato = readCbModel([projDir 'Data/comparison_models/' 'curatedPotatoGEM.xml']);
Sobliquus = readCbModel([projDir 'Data/comparison_models/' 'iAR632.xml']);
Cohadii = readCbModel([projDir 'Data/comparison_models/' 'iCO1515_mixo.xml']);
aracore = readCbModel([projDir 'Data/comparison_models/'...
    'AraCore_v2_1.wKEGG.xml']);
sweetlove = readCbModel([projDir 'Data/comparison_models/' ...
    'CAM_diel_sweetlove_valid.xml']);
soy = readCbModel([projDir 'Data/comparison_models/' 'soybean_wComps.xml']);
pcm = readCbModel([projDir 'Data/pcm/' 'pcm.v' num2str(pcmImportVersion) '.xml']);

nModels = 13;
%% extract chloroplast networks for each model

data.model_name = cell(nModels, 1);
data.nChlRxns = zeros(nModels, 1);
data.totModelRxns = zeros(nModels, 1);
data.nTransp = zeros(nModels, 1);
data.nIntraChlTransp = zeros(nModels, 1);
data.nMetsStroma = zeros(nModels, 1);
data.nMetsThylLum = zeros(nModels, 1);
data.nMetsChl = zeros(nModels, 1);
data.nGenes = zeros(nModels, 1);
data.nSubsystems = zeros(nModels, 1);
data.nChloroplastComps = zeros(nModels, 1);
data.chlComps = cell(nModels, 1);
data.nPanGenes = zeros(nModels, 1);

modelIx = 1;

% aracore
chlComps = {'h', 'l'};
metComps = cell(length(aracore.mets), 1);
for i = 1:length(aracore.mets)
    res = regexp(aracore.mets{i}, '\[(.*?)\]', 'tokens');
    if ~isempty(res)
        metComps(i) = res{end};
    end
end

data.model_name{modelIx} = 'aracore';
data = gather_data(data, aracore, chlComps, 'h', 'l', ...
    metComps, modelIx, ',\s');

modelIx = modelIx + 1;

% aragem
chlComps = {'p'};
metComps = cell(length(aragem.mets), 1);
for i = 1:length(aragem.mets)
    res = regexp(aragem.mets{i}, '\[(.*?)\]', 'tokens');
    if ~isempty(res)
        metComps(i) = res{end};
    end
end

data.model_name{modelIx} = 'aragem';
data = gather_data(data, aragem, chlComps, 'p', 'non_existent', ...
    metComps, modelIx);

modelIx = modelIx + 1;

% sweetlove (diel model, use both parts)
chlComps = {'p1', 'p2', 'l1', 'l2'};
metComps = cell(length(sweetlove.mets), 1);
for i = 1:length(sweetlove.mets)
    res = regexp(sweetlove.mets{i}, '\[(.*?)\]', 'tokens');
    if ~isempty(res)
        metComps(i) = res{end};
    end
end

data.model_name{modelIx} = 'sweetlove';
data = gather_data(data, sweetlove, chlComps, {'p1', 'p2'}, {'l1', 'l2'}, ...
    metComps, modelIx);

modelIx = modelIx + 1;

% maize
chlComps = {'d0'};
metComps = cell(length(maize.mets), 1);
for i = 1:length(maize.mets)
    res = regexp(maize.mets{i}, '\[(.*?)\]', 'tokens');
    if ~isempty(res)
        metComps(i) = res{end};
    end
end

data.model_name{modelIx} = 'maize';
data = gather_data(data, maize, chlComps, {'d0'}, 'non-existent', ...
    metComps, modelIx);

modelIx = modelIx + 1;

% rice
chlComps = {'s', 'u'};
metComps = cell(length(rice.mets), 1);
for i = 1:length(rice.mets)
    res = regexp(rice.mets{i}, '\[(.*?)\]', 'tokens');
    if ~isempty(res)
        metComps(i) = res{end};
    end
end

data.model_name{modelIx} = 'rice';
data = gather_data(data, rice, chlComps, 's', 'u', metComps, modelIx);

modelIx = modelIx + 1;

% mintzOron
chlComps = {'Plastid'};
metComps = cell(length(mintzOron.mets), 1);
for i = 1:length(mintzOron.mets)
    res = regexp(mintzOron.mets{i}, '\[(.*?)\]', 'tokens');
    if ~isempty(res)
        metComps(i) = res{end};
    end
end

data.model_name{modelIx} = 'mintzOron';
data = gather_data(data, mintzOron, chlComps, 'Plastid', 'non-existent', ...
    metComps, modelIx, '(?<!\s)/(?!\s)'); % match '/' without spaces on either side

modelIx = modelIx + 1;

% panAlgae
chlComps = {'h', 's', 'u'};
metComps = cell(length(panAlgae.mets), 1);
for i = 1:length(panAlgae.mets)
    res = regexp(panAlgae.mets{i}, '\[(.*?)\]', 'tokens');
    if ~isempty(res)
        metComps(i) = res{end};
    end
end

data.model_name{modelIx} = 'panAlgae';
data = gather_data(data, panAlgae, chlComps, 'h', 'u', ...
    metComps, modelIx);

modelIx = modelIx + 1;

% poplar
chlComps = {'p', 't', 'i'};
metComps = cell(length(poplar.mets), 1);
for i = 1:length(poplar.mets)
    res = regexp(poplar.mets{i}, '\[(.*?)\]', 'tokens');
    if ~isempty(res)
        metComps(i) = res{end};
    end
end

data.model_name{modelIx} = 'poplar';
data = gather_data(data, poplar, chlComps, 'p', 't', ...
    metComps, modelIx);

modelIx = modelIx + 1;

% soy
chlComps = {'p'};
metComps = cell(length(soy.mets), 1);
for i = 1:length(soy.mets)
    res = regexp(soy.mets{i}, '\[(.*?)\]', 'tokens');
    if ~isempty(res)
        metComps(i) = res{end};
    end
end

data.model_name{modelIx} = 'soy';
data = gather_data(data, soy, chlComps, 'p', 'non-existent', ...
    metComps, modelIx);

modelIx = modelIx + 1;

% potato
% has no compNames field (it does but is equal to comps), so I had to make
% some educated guesses based on compartment names of the PLM by Sandra 
% Correa, who also participated in this paper
chlComps = {'h', 'hm', 'ohm', 'pg', 'l'};
metComps = cell(length(potato.mets), 1);
for i = 1:length(potato.mets)
    res = regexp(potato.mets{i}, '\[(.*?)\]', 'tokens');
    if ~isempty(res)
        metComps(i) = res{end};
    end
end

data.model_name{modelIx} = 'potato';
data = gather_data(data, potato, chlComps, 'h', 'l', ...
    metComps, modelIx);

modelIx = modelIx + 1;

% Sobliquus
chlComps = {'d0'};
metComps = cell(length(Sobliquus.mets), 1);
for i = 1:length(Sobliquus.mets)
    res = regexp(Sobliquus.mets{i}, '\[(.*?)\]', 'tokens');
    if ~isempty(res)
        metComps(i) = res{end};
    end
end

data.model_name{modelIx} = 'Sobliquus';
data = gather_data(data, Sobliquus, chlComps, 'd0', 'non-existent', ...
    metComps, modelIx);

modelIx = modelIx + 1;

% Cohadii
chlComps = {'ch', 'th'};
metComps = cell(length(Cohadii.mets), 1);
for i = 1:length(Cohadii.mets)
    res = regexp(Cohadii.mets{i}, '\[(.*?)\]', 'tokens');
    if ~isempty(res)
        metComps(i) = res{end};
    end
end

data.model_name{modelIx} = 'Cohadii';
data = gather_data(data, Cohadii, chlComps, 'ch', 'th', ...
    metComps, modelIx);

modelIx = modelIx + 1;

% pcm
chlComps = {'h', 'hm', 'ohm', 'pg', 'l'};
metComps = cell(length(pcm.mets), 1);
for i = 1:length(pcm.mets)
    res = regexp(pcm.mets{i}, '\[(.*?)\]', 'tokens');
    if ~isempty(res)
        metComps(i) = res{end};
    end
end

data.model_name{modelIx} = 'pcm';
data = gather_data(data, pcm, chlComps, 'h', 'l', ...
    metComps, modelIx, '(?<!\s),(?!\s)');

% get number of union genes in pcm
disp(['Number of union genes in pcm: ' num2str(data.nGenes(modelIx))]);

% get number of genes in any species-specfic version of the pcm, and of
% Arabidopsis specifically
modelsDir = [projDir 'Data/pcm/species/'];
currPattern = ['*.pcm.v', num2str(pcmImportVersion), '.mat'];
fileList = dir(fullfile(modelsDir, currPattern));
orgs = cell(length(fileList), 1);

unionGenes = pcm.genes;
pcmGenes = unionGenes;
for i = 1:length(fileList)
    % Remove the suffix to get the base name
    orgs{i} = strtok(fileList(i).name, '.');
    newerVersion = strrep(fileList(i).name, ...
            num2str(pcmImportVersion), num2str(pcmImportVersion + 1));
    if isfile([modelsDir newerVersion])
        load([modelsDir newerVersion]);
    else
        load([modelsDir fileList(i).name]);
    end
    pcmGenes = unique([pcmGenes; setdiff(model.genes, unionGenes)]);

    if contains(orgs{i}, 'Arabidopsis')
        disp(['Number of genes in Ara pcm: ' num2str(length(model.genes))]);
    end
end
data.nPanGenes(modelIx) = length(unique(pcmGenes));
disp(['Number of genes of any species in pcm: ' num2str(data.nPanGenes(modelIx))]);

table = struct2table(data);
writetable(table, [projDir 'Data/analysis/modelComparison/modelComparison.csv']);

%% Functions
function [chlRxns, chlTransp, intraChlTransp] = extractRxns(...
    model, chlComps, metComps)
    chlRxnsIx = zeros(1, length(model.rxns));
    chlTranspIx = zeros(1, length(model.rxns));
    intraChlTranspIx = zeros(1, length(model.rxns));
    
    for i = 1:length(model.rxns)
        
        % Get the compartments of each metabolite involved
        currMetComps = unique(metComps(model.S(:, i) ~= 0));
        
        % Check if all metabolites are in a single specified compartment
        if all(ismember(currMetComps, chlComps))
            % If all in the same one, class as chlRxns
            if isscalar(unique(currMetComps))
                chlRxnsIx(i) = 1;
            else % All metabolites are in chlcomps but in different compartments
                intraChlTranspIx(i) = 1;
            end
        elseif any(ismember(currMetComps, chlComps))
            % Metabolites in some but not all in chlcomps -> transporter
            chlTranspIx(i) = 1;
        end
    end

    chlRxns = model.rxns(logical(chlRxnsIx));
    chlTransp = model.rxns(logical(chlTranspIx));
    intraChlTransp = model.rxns(logical(intraChlTranspIx));
end

% extract nMets of stroma, thylakoid/lumen, and plastid total
function [stromaMets, thylMets, totalMets] = ...
    extractMets(model, metComps, stromaComp, thylakoidComp, chlComps)

    % initialize indices vectors
    stromaMetsIx = zeros(1, length(metComps));
    thylMetsIx = zeros(1, length(metComps));
    totalMetsIx = zeros(1, length(metComps));

    % for all metabolites, classify their compartment
    for i = 1:length(metComps)
        comp = metComps{i};
        if any(strcmp(comp, stromaComp))
            stromaMetsIx(i) = 1;
            totalMetsIx(i) = 1;
        elseif any(strcmp(comp, thylakoidComp))
            thylMetsIx(i) = 1;
            totalMetsIx(i) = 1;
        elseif any(strcmp(chlComps, comp))
            totalMetsIx(i) = 1;
        end
    end

    % assemble the output vectors
    stromaMets = model.mets(logical(stromaMetsIx));
    thylMets = model.mets(logical(thylMetsIx));
    totalMets = model.mets(logical(totalMetsIx));
end

function data = gather_data(data, model, chlComps, stromaComp, ...
    thylComp, metComps, modelIx, subsRegexpSplit)

    data.totModelRxns(modelIx) = length(model.rxns);
    
    % extract nRxns of chl, chl/outside transporters, intrachloro transporters
    [chlRxns, chlTransp, intraChlTransp] = ...
        extractRxns(model, chlComps, metComps);
    allChlRxns = [chlRxns; chlTransp; intraChlTransp];
    data.nChlRxns(modelIx) = length(chlRxns);
    data.nTransp(modelIx) = length(chlTransp);
    data.nIntraChlTransp(modelIx) = length(intraChlTransp);
    
    % extract nMets of stroma, thylakoid/lumen, and plastid total
    [stromaMets, thylMets, totalMets] = ...
        extractMets(model, metComps, stromaComp, thylComp, chlComps);
    data.nMetsStroma(modelIx) = length(stromaMets);
    data.nMetsThylLum(modelIx) = length(thylMets);
    data.nMetsChl(modelIx) = length(totalMets);
    
    % nGenes (model.genes entails all genes of the model, not only those
    % acting in the chloroplast)
    if isfield(model, 'rules')
        genesIx = zeros(1, length(model.genes));
        for j = 1:length(allChlRxns)
            res = regexp(model.rules(strcmp(model.rxns, allChlRxns{j})), ...
                                     'x\((.*?)\)', 'tokens');
            if ~isempty(res)
                for i = 1:length(res{1})
                    genesIx(str2double(char(res{1}{i}))) = 1;
                end
            end
        end
    else
        disp(modelIx);
    end
    data.nGenes(modelIx) = sum(genesIx);
    
    % nSubsystems
    if isfield(model, 'subSystems')
        if iscell(model.subSystems{1})
            rxnIx = find(ismember(model.rxns, chlRxns) | ...
                    ismember(model.rxns, intraChlTransp));
            chlSubsystems = {};
            for rxnIxIx = 1:length(rxnIx)
                currSubs = model.subSystems{rxnIx(rxnIxIx)};
                for subsIx = 1:length(currSubs)
                    chlSubsystems = unique([chlSubsystems, ...
                                            currSubs(subsIx)]);
                end
            end
        else
            chlSubsystems = getUniqueSubs(model.subSystems(...
                            ismember(model.rxns, chlRxns) | ...
                            ismember(model.rxns, intraChlTransp)), ...
                            subsRegexpSplit);
        end
        chlSubsystems = unique(cellfun(@lower, ...
            chlSubsystems(~cellfun(@isempty, chlSubsystems)).', ...
            'UniformOutput',false));
        data.nSubsystems(modelIx) = length(chlSubsystems);
    else
        data.nSubsystems(modelIx) = 0;
    end
    
    % nChloroplastComps
    data.nChloroplastComps(modelIx) = length(chlComps);
    
    % chlComps
    data.chlComps(modelIx) = {chlComps};
end

function uniqueSubs = getUniqueSubs(subs, regexpPattern)
    % input:  subs: cell array containing strings
    %         regexp: regexp that matches sequences where individual strings
    %                are split
    % output: unique cell array containing each sub-element of subs after
    %         subs was split

    allSubs = {};

    for i = 1:length(subs)
        if ~isempty(subs{i})
            % Use regexp to split the current string with the given pattern
            splitElements = regexp(subs{i}, regexpPattern, 'split');
            allSubs = [allSubs, splitElements]; %#ok<AGROW>
        end
    end

    uniqueSubs = unique(allSubs);
end

function KEGGs = extractKeggPanAlgae(rxnNotes)
    % extract the KEGG IDs from the rxnNotes field of the panAlgae model

    KEGGs = cell(length(rxnNotes), 1);
    for i = 1:length(rxnNotes)
        currNote = rxnNotes{i};
        if isempty(currNote)
            continue
        end

        % the rxnNotes can contain multiple lines; KEGG IDs would be on
        % their own and be preceeded by 'References:'
        lines = strsplit(currNote, newline);
        currKeggs = char.empty;
        for j = 1:length(lines)
            if ~contains(lines{j}, 'References:')
                continue
            end
            currKeggs = strrep(lines{j}, 'References:', '');
        end

        % few lines have dois as their reference
        if ~isempty(currKeggs) && ~contains(currKeggs, 'doi')
            KEGGs{i} = strrep(currKeggs, ';', ',');
        end
    end
end