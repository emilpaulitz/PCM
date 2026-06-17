%% Generate a organism-specific chloroplast model(s)
% The script generates organism-specific PCM versions for all organisms for
% which BLAST results are found in blastDir (see below)
clearvars -except gurobiAvailable projDir; clc;

%% ADAPT THESE:
% when mapping complex gprs, do we require all subunits in unionPCM to have
% an ortholog in the new organism to consider the gpr valid?
requireAllSubunits = false;
% if empty, will use all organisms with BLAST results in blastDir
listOfOrgs = {'Bay-0', 'Bla-1'};
% folder relative to projDir with BLAST result files named ${ID}.out.tsv
% (see Code/accSpecPCM/blast.sh for example commands)
blastDir = 'Data/analysis/accSpecPCM/03_BLAST_files/';
% folder relative to projDir in which the resulting models will be written
modelOutPath = 'Data/analysis/accSpecPCM/models/';
% the format(s) in which the model should be written out
writeMat = true;
writeXML = true;
% version of union pcm to use
importVersion = 1;

%% Start processing
if ~writeMat && ~writeXML
    error(['Specfiy one output format to be written, otherwise the ' ...
        'script does nothing!']) %#ok<UNRCH>
end
% read media definitions
media = readtable([projDir 'Data/analysis/simulationMedia.tsv'], ...
    'FileType', 'text', 'Delimiter', '\t');
media.hetero = media.hetero_str;
media = removevars(media, ["hetero_suc", "hetero_str"]);

% load union PCM model
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

% get demand production rate under phototrophic conditions
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

% gather list of all accession names
if isempty(listOfOrgs)
    fileList = dir(fullfile(blastDir, '*.out.tsv'));
    listOfOrgs = cell(length(fileList), 1);
    for i = 1:length(fileList)
        % Remove the suffix to get the base name
        [~, baseName, ~] = fileparts(fileList(i).name);
        listOfOrgs{i} = baseName(1:end-length('.out'));
    end
end

%% make org specific: replace gprs, find missing reactions
for orgIx = 1:length(listOfOrgs)
    org = listOfOrgs{orgIx};
    
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    disp(['Working on ' org])
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

    % make org specific: replace gprs, find missing reactions
    [orgPcm, missingInOrg] = makeOrgSpecific(pcm, ...
        [projDir, blastDir, org, '.out.tsv'], requireAllSubunits);
    solInit = solveLP(orgPcm);

    % set up variables for carveMe procedure
    M = 1000;  % upper/lower bound
    eps = 1e-6;  % min flux
    % min expected growth; it would be unrealistic to have much less
    % production of biomass components than in another species/accessions
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

    % solve; if gurobi is used and the minimum growth is not reached,
    % increase feasibility and try again
    while ~solved && ~broken

        [y_res, v_res] = carvePCM(orgPcm, M, eps, ...
            min_gr_photo, min_gr_hetero, ... 
            rxn_scores, media, gurobiAvailable, true, FeasibilityTol);
        
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
    
    if isfield(carved, 'A')
        carved = rmfield(carved, 'A');
    end
    if writeMat
        writeCbModel(carved, 'fileName', [modelPath, 'species/', org, ...
                                '.pcm.v' num2str(importVersion) '.mat']);
    end
    if writeXML
        writeCbModel(carved, 'fileName', [modelPath, 'species/', org, ...
                                '.pcm.v' num2str(importVersion) '.xml']);
    end
end
