%% Generate species-specific chloroplast models
% Use as genetic evidence both the EC-gene associations and BLAST results
clearvars -except gurobiAvailable projDir; clc;

% when mapping complex gprs, do we require all Arabidopsis subunits to have
% an ortholog in the new organism to consider the gpr valid?
requireAllSubunits = false;
importVersion = 1;

% read media definitions
media = readtable([projDir 'Data/analysis/simulationMedia.tsv'], ...
    'FileType', 'text', 'Delimiter', '\t');
media.hetero = media.hetero_str;
media = removevars(media, ["hetero_suc", "hetero_str"]);

% load model
modelPath = strcat(projDir, 'Data/pcm/');

% get demand production rate for heterotrophic conditions
load(strcat(modelPath, 'pcm.v', num2str(importVersion), '.mat'), 'model');
pcm = model;
for rxnIx = 1:length(media.rxn)
    currRxnIx = strcmp(media.rxn{rxnIx}, model.rxns);
    if any(currRxnIx)
        pcm.ub(currRxnIx) = media.hetero(rxnIx);
    end
end
gr_het = solveLP(pcm).f;

% get demand production rate phototrophic conditions
load(strcat(modelPath, 'pcm.v', num2str(importVersion), '.mat'), 'model');
pcm = model;
for rxnIx = 1:length(media.rxn)
    currRxnIx = strcmp(media.rxn{rxnIx}, model.rxns);
    if any(currRxnIx)
        pcm.ub(currRxnIx) = media.photo(rxnIx);
    end
end
gr_photo = solveLP(pcm).f;

% load the original again
load(strcat(modelPath, 'pcm.v', num2str(importVersion), '.mat'), 'model');
pcm = model;

% gather list of all orgs
fileList = dir(fullfile([projDir 'Data/sequences/'], '*.fasta'));
orgs = cell(length(fileList), 1);
for i = 1:length(fileList)
    % Remove the suffix to get the base name
    [~, baseName, ~] = fileparts(fileList(i).name);
    orgs{i} = baseName;
end

%% make org specific: replace gprs, find missing reactions
for orgIx = length(orgs):-1:1
    org = orgs{orgIx};

    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    disp(['Working on ' org])
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

    % make org specific: replace gprs, find missing reactions
    [orgPcm, missingInOrg, orgRxnGenes] = makeSpeciesSpecificECAssoc(...
        pcm, org, requireAllSubunits, projDir);
    solInit = solveLP(orgPcm);

    % set up variables for carveMe procedure
    M = 1000;  % upper/lower bound
    eps = 1e-6;  % min flux
    % min expected growth; it would be unrealistic to have much less
    % production of biomass components than in another species
    min_gr_photo = 0.8 * gr_photo;
    min_gr_hetero = 0.8 * gr_het;

    rxn_scores = ones(length(orgPcm.rxns), 1);
    rxn_scores(logical(missingInOrg)) = -1;
    rxn_scores(strcmpi(pcm.grRules, 'unknown') | ...
        strcmpi(pcm.grRules, 'brenda evidence')) = 0;
    rxn_scores(contains(pcm.subSystems, 'pseudo') | ...
               contains(pcm.subSystems, 'import') | ...
               contains(pcm.subSystems, 'exchange') | ...
               contains(pcm.subSystems, 'export')) = 0;
    rxn_scores(contains(pcm.subSystems, 'transport')) = 1;
    solved = false;
    broken = false;
    FeasibilityTol = 1e-6;

    % solve; if the minimum growth is not reached, increase feasibility and
    % try again
    while ~solved && ~broken

        [y_res, v_res] = carvePCM(orgPcm, M, eps, min_gr_photo, ...
            min_gr_hetero, rxn_scores, media, gurobiAvailable, true, ...
            FeasibilityTol);
        
        % Carve model
        y = y_res(1:length((pcm.rxns)));
        rxnsToRemove = pcm.rxns(y < 1e-5);
        carved = removeRxns(orgPcm, rxnsToRemove);
        sol = solveLP(carved);
        disp(['Initial demand rxn flux: ', num2str(solInit.f)]);
        disp(['Initial demand rxn flux under phototropic conditions: ', ...
            num2str(gr_photo)]);
        disp(['Initial demand rxn flux under heterotropic conditions: ',...
            num2str(gr_het)]);

        % Check biomass production under different conditions
        oriUb = carved.ub;
        disp(['Carved model demand rxn flux: ', num2str(sol.f)]);
        for rxnIx = 1:length(media.rxn)
            currRxnIx = strcmp(media.rxn{rxnIx}, carved.rxns);
            if any(currRxnIx)
                carved.ub(currRxnIx) = media.photo(rxnIx);
            end
        end
        carvedPhotoSol = solveLP(carved);
        disp(['Carved model demand rxn under phototropic conditions: ', ...
            num2str(carvedPhotoSol.f)]);

        for rxnIx = 1:length(media.rxn)
            currRxnIx = strcmp(media.rxn{rxnIx}, carved.rxns);
            if any(currRxnIx)
                carved.ub(currRxnIx) = media.hetero(rxnIx);
            end
        end
        carvedHeteroSol = solveLP(carved);
        disp(['Carved model demand rxn under heterotropic conditions: ', ...
            num2str(carvedHeteroSol.f)]);
        carved.ub = oriUb;

        % handle the different cases if the biomass constraint is not
        % fulfilled anymore
        if carvedPhotoSol.f >= min_gr_photo - 1e-8 || ...
                carvedHeteroSol.f >= min_gr_hetero - 1e-8
            solved = true;
        elseif gurobiAvailable
            if FeasibilityTol >= 1e-8
                FeasibilityTol = 1e-9;
                eps = 1e-8;
                disp(['Solution not found, reducing FeasibilityTol to ' ...
                    num2str(FeasibilityTol) ' and eps to ' num2str(eps)]);
            else
                disp(['Solution not found but FeasibilityTol cannot be ' ...
                    'reduced further']);
                broken = true;
            end
        else 
            disp('No solution found that can produce biomass!');
            broken = true;
        end
    end
    if broken
        disp(['Did not find a solution for ' org])
        continue
    end

    % if we found a working solution, produce and write out the model
    carved = updateFromGrRules(carved);
    
    % write out model
    if isfield(carved, 'A')
        carved = rmfield(carved, 'A');
    end
    writeCbModel(carved, 'fileName', ...
        outModelPath(modelPath, org, importVersion));
end

function outPath = outModelPath(modelPath, org, importVersion)
    outPath = [modelPath, 'species/', org, '.pcm.v'...
               num2str(importVersion) '.mat'];
end


%% main function
% like makeOrgSpec.m, but also integrates evidence from the original
% prediction of EC numbers
function [orgPcm, missingInOrg, orgRxnGenes] = makeSpeciesSpecificECAssoc(pcm, org, ...
            requireAllSubunits, projPath)

    % Thresholds for BLAST results
    identityThreshold = 80;
    seqLengthThreshold = 2.0;
    hitLengthThreshold = 0.5;

    % read files
    orgECFile = [projPath 'Data/EC_predictions/PredBy2Tool/' org '.tsv'];
    orgBLASTres = [projPath 'Data/BLAST_results_37_species/' org ...
        '.out.tsv'];
    
    % Load EC prediction results
    orgECData = readtable(orgECFile, 'FileType', 'text', 'Delimiter', '\t');
    orgGenes = orgECData.Gene;
    orgEcNumbers = orgECData.EC;

    % Load BLAST results
    blastResults = readtable(orgBLASTres, 'FileType', 'text', ...
        'Delimiter', '\t');
    
    tableOrgGenes = blastResults.Var1;
    tableUnionGenes = blastResults.Var2;
    tablePids = blastResults.Var3;
    tableQlen = blastResults.Var7;
    tableSlen = blastResults.Var8;
    tableHitlen = blastResults.Var9;

    % Use EC-associations of the species to associate genes to reactions in
    % the model
    orgRxnGenes = cell(length(pcm.rxns), 1);
    for i = 1:length(pcm.rxns)
        % Initialize an empty cell array to store associated genes
        associatedGenes = {}; 
        
        pcmECList = strsplit(pcm.rxnEC{i}, ',');
        for j = 1:length(pcmECList)
            currPcmEC = pcmECList{j};
    
            % Split each EC-number from the organism file and check for matches
            for k = 1:length(orgEcNumbers)
                ecNumberList = strsplit(orgEcNumbers{k}, '; ');
                
                % Check if currPcmEC is in the list of EC numbers
                if any(strcmp(ecNumberList, currPcmEC))
                    % Add the corresponding gene to the list of associated genes
                    associatedGenes = [associatedGenes; ...
                                       strrep(orgGenes{k}, '|', '_')];
                end
            end
        end
    
        orgRxnGenes{i} = unique(associatedGenes);
    end
    
    % Gather species-genes from BLAST results
    
    % Obtain a mask that filters relevant hits
    seqLenRatios = tableQlen ./ tableSlen;
    relevantBlastHit = (tablePids >= identityThreshold) & ...
        (seqLenRatios <= seqLengthThreshold) & ...
        (seqLenRatios >= (1/seqLengthThreshold)) & ...
        (tableHitlen ./ min(tableSlen, tableQlen) >= hitLengthThreshold);

    % copy the union model to the organism model
    orgPcm = pcm;
    
    % an array to track reactions missing genes in the species
    missingInOrg = false(length(orgPcm.rxns), 1);
    
    for i = 1:length(orgPcm.rxns)
        currentUnionGpr = pcm.grRules{i};
    
        % Check for reactions without any genes (since we anyway only
        % associate species genes to union model genes in this step) and
        % non-enzymatic rxns
        if isempty(pcm.grRules{i}) || ...
                any(strcmpi(currentUnionGpr, {'pseudo', 'spontaneous', 'unknown'}))
            continue;
        end
    
        unionGenes = pcm.genes(logical(pcm.rxnGeneMat(i, :)));
        currRxnOrgGenes = {};
        
        % Use BLAST results to map genes
        for j = 1:length(unionGenes)
            % Find rows in BLAST results corresponding to the union gene
            idxBlast = strcmpi(tableUnionGenes, unionGenes{j}) & relevantBlastHit;
            
            % Map organism genes based on BLAST results
            mappedGenesTmp = tableOrgGenes(idxBlast);
            if ~isempty(mappedGenesTmp)
                currRxnOrgGenes = [currRxnOrgGenes; ...
                               strrep(mappedGenesTmp, '|', '_')];
            end
        end
    
        orgRxnGenes{i} = unique(vertcat(orgRxnGenes{i}, currRxnOrgGenes));
    end
    
    % Build species-specific GPR-rules, based on found associations

    % Update grRules with new mappings
    for i = 1:length(orgPcm.rxns)

        if any(strcmpi(currentUnionGpr, ...
                {'pseudo', 'spontaneous', 'unknown'}))
            continue;
        end

        currentUnionGpr = pcm.grRules{i};
    
        if isempty(orgRxnGenes{i})

            % Mark reaction in the missingInOrg array if no mapping found
            missingInOrg(i) = true;

            orgPcm.grRules{i} = char.empty;
    
        elseif ~contains(currentUnionGpr, ' AND ', 'IgnoreCase', true)
            % AND is more complicated 
            orgPcm.grRules{i} = strjoin(unique(orgRxnGenes{i}), ' Or ');
    
        else
            % in case of AND, a one-to-one mapping of union model genes to species 
            % genes is required. Thus, any genes found via EC-gene associations
            % that are not found by BLAST, can not be used automatically and 
            % require manual curation
    
            % these two arrays will hold the union model genes associated to
            % the reaction, and the respective orthologs of the current species
            unionGenes = pcm.genes(logical(pcm.rxnGeneMat(i, :)));
            replacements = cell(length(unionGenes), 1);
            
            for j = 1:length(unionGenes)
                % Find rows in BLAST results corresponding to the Arabidopsis gene
    
                idxBlast = strcmpi(tableUnionGenes, unionGenes{j}) & relevantBlastHit;
                
                % Map organism genes based on BLAST results
                mappedGenesTmp = tableOrgGenes(idxBlast);
                if ~isempty(mappedGenesTmp)
                    orgGeneMapping = strjoin(unique(mappedGenesTmp), ' Or ');
                    orgGeneMapping = strrep(orgGeneMapping, '|', '_');
                    replacements{j} = ['( ' orgGeneMapping ' )'];
                end
            end
                
            % replace genes 1 to 1
            % first replace by placeholder to avoid problems if union gene
            % names occur in new gene names
            orgGprOneToOne = pcm.grRules{i};
            placeholder = 'xx$$%%';
            for j = 1:length(unionGenes)
                orgGprOneToOne = strrep(orgGprOneToOne, unionGenes{j}, ...
                                        [placeholder '[' num2str(j) ']']);
            end
    
            if requireAllSubunits

                % then replace placeholder by new gene(s) or empty string
                for j = 1:length(unionGenes)
                    if isempty(replacements{j})
                        replacements{j} = char.empty;
                    end
                    orgGprOneToOne = strrep(orgGprOneToOne, ...
                        [placeholder '[' num2str(j) ']'], replacements{j});
    
                end
    
                % clean up the gpr, leaving it empty if a subunit is
                % missing
                orgGprOneToOne = format_gpr(orgGprOneToOne);

                if isempty(orgGprOneToOne)
                    missingInOrg(i) = true;
                end
                % Create a new GPR rule using the newly mapped genes, leaving 
                % a hole in the gpr if no mappings exist
                orgPcm.grRules{i} = orgGprOneToOne;
            else
                % replace placeholder by new gene(s), leaving union genes 
                % if no mappings exist
                for j = 1:length(unionGenes)
                    if isempty(replacements{j})
                        orgGprOneToOne = strrep(orgGprOneToOne, ...
                            [placeholder '[' num2str(j) ']'], ...
                            unionGenes{j});
                    else
                        orgGprOneToOne = strrep(orgGprOneToOne, ...
                            [placeholder '[' num2str(j) ']'], ...
                            replacements{j});
                    end
                end
        
                orgPcm.grRules{i} = orgGprOneToOne;
            end
        end
    end
end

%% Helper functions
% these following functions are translated from my python functions with
% ChatGPT, but seem to work fine
% format_gpr handles all cases, removing any isoenzymes that are not
% functional because a subunit is missing. That can mean the entire string
% is empty if a subunit of the only, or all, isozymes is missing
function formatted_gpr = format_gpr(gpr)
    formatted_gpr = simplify_gpr(gpr);
    if contains(lower(gpr), 'and')
        formatted_gpr = format_and_gpr(formatted_gpr);
    end
end

function tokens = split_gpr(gpr)
    tokens = {};
    words = strsplit(gpr, ' ');  % Split the GPR by spaces
    
    for i = 1:length(words)
        word = words{i};
        if strcmpi(word, 'or') || strcmpi(word, 'and')
            tokens{end+1} = word; %#ok<*AGROW> % Add logical operator directly
        else
            while contains(word, '(')
                word = extractAfter(word, '(');
                tokens{end+1} = '(';  % Include opening bracket
            end
            
            num_close_brackets = 0;
            while contains(word, ')')
                word = extractBefore(word, ')');
                num_close_brackets = num_close_brackets + 1;
            end
            
            if ~isempty(word)
                tokens{end+1} = word;  % Append non-bracket word token
            end
            
            % Append closing brackets accordingly
            tokens(end+1:end+num_close_brackets) = repmat({')'}, 1, num_close_brackets);
        end
    end
end

function simplified_gpr = simplify_gpr(curr_gpr)
    curr_gpr = strtrim(curr_gpr);  % Trim leading and trailing whitespace
    changed = true;
    
    while changed
        changed = false;
        curr_gpr = strtrim(curr_gpr);
        
        % Handle multiple logical cleaning and orphan bracket issues
        if startsWith(lower(curr_gpr), 'or ')
            curr_gpr = curr_gpr(4:end);
            changed = true;
        end
        if endsWith(lower(curr_gpr), ' or')
            curr_gpr = curr_gpr(1:end-3);
            changed = true;
        end
        if contains(curr_gpr, 'or or')
            curr_gpr = strrep(curr_gpr, 'or or', 'or');
            changed = true;
        end
        if contains(curr_gpr, 'OR OR')
            curr_gpr = strrep(curr_gpr, 'OR OR', 'OR');
            changed = true;
        end
        if contains(curr_gpr, '  ')
            curr_gpr = strrep(curr_gpr, '  ', ' ');
            changed = true;
        end
        if contains(curr_gpr, '( or')
            curr_gpr = strrep(curr_gpr, '( or', '(');
            changed = true;
        end
        if contains(curr_gpr, '( OR')
            curr_gpr = strrep(curr_gpr, '( OR', '(');
            changed = true;
        end
        if contains(curr_gpr, 'OR )')
            curr_gpr = strrep(curr_gpr, 'OR )', ')');
            changed = true;
        end
        if contains(curr_gpr, '( )')
            curr_gpr = strrep(curr_gpr, '( )', '');
            changed = true;
        end
    end
    
    simplified_gpr = curr_gpr;
end

function formatted_gpr = format_and_gpr(gpr)

    gpr = strtrim(gpr);

    % gpr is not functional
    if contains(lower(gpr), 'and and') ||...
            endsWith(lower(gpr), 'and') || startsWith(lower(gpr), 'and')
        formatted_gpr = '';
    else
        % refer complicated cases
        formatted_gpr = format_and_gpr_worker(gpr);
    end
end

function formatted_gpr = format_and_gpr_worker(s)
    prev_eles = {};
    curr_eles = {};
    curr_depth = 0;
    rem_at_depth = false;
    words = split_gpr(s);  % Use the split_gpr defined above

    prev = '';
    for i = 1:length(words)
        word = words{i};
        
        if rem_at_depth
            if strcmp(word, ')')
                curr_depth = curr_depth - 1;
            elseif strcmp(word, '(')
                curr_depth = curr_depth + 1;
            end
            
            if curr_depth >= rem_at_depth
                prev = '';
                continue;
            else
                rem_at_depth = false;
            end
        end
        
        if strcmp(word, '(')
            curr_depth = curr_depth + 1;
            prev_eles{end+1} = curr_eles; %#ok<*AGROW>
            curr_eles = {};
        end
        
        if strcmpi(word, 'and') && strcmp(prev, '(')
            if ~rem_at_depth
                rem_at_depth = curr_depth;
            end
            continue;
        end
        
        if strcmp(word, ')')
            if strcmpi(prev, 'and')
                % Discard curr_eles because the current element is not functional
                curr_eles = prev_eles{end};
                prev_eles(end) = [];
                
                % Do not append the current ')' to curr_eles
                prev = '';
                continue;
            else
                curr_eles = [prev_eles{end}, curr_eles];
                prev_eles(end) = [];
            end
        end
        
        curr_eles{end+1} = word;
        prev = word;
    end
    
    formatted_gpr = simplify_gpr(strjoin(curr_eles, ' '));  % Use simplify_gpr defined above
end