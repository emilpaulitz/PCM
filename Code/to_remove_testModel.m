clearvars -except gurobiAvailable projDir; clc;

modelPath = strcat(projDir, 'Data/pcm/');
importVersion = 1;
load([modelPath 'pcm.v' num2str(importVersion) '.mat'], 'model');

%% 
% Update metComps array
for i = 1:length(model.mets)
    model.metComps(i) = find(strcmp(model.comps, extractComp(model.mets{i})));
end

model = buildRxnEquations(model);

%% How many transporters do we have?
% Get transporters to the cytosol for inspection
cytIdx = find(strcmp(model.comps, 'c'));
cytTranspIdx = [];

lastPLMIdx = find(strcmp(model.rxns, 'Tr_tRNA'));
numAfterPLM = 0;
for i = 1:length(model.rxns)
    % Get metabolite indices associated with the current reaction
    met_indices = find(model.S(:, i) ~= 0);

    % Check compartments for these metabolites
    comps = unique(model.metComps(met_indices));

    % If 'c' compartment is present along with other compartments
    if ismember(cytIdx, comps) && length(comps) > 1 && ...
        ~startsWith(model.rxns(i), 'Exch_') && ~strcmp(model.rxns(i), 'precursorPool')
        cytTranspIdx = [cytTranspIdx; i];

        if i > lastPLMIdx
            numAfterPLM = numAfterPLM + 1;
        end
    end
end

disp('Number of cytosol transporters:')
disp(length(cytTranspIdx))

disp('Number of cytosol transporters from my curation:')
disp(numAfterPLM)

% If there is something in the notes field, it comes from the curation by
% Robin
% Otherwise, it is from the PLM. Two reactions without notes are also from
% Robin
rxnIdxEmptyNotes = cytTranspIdx(cellfun(@isempty, model.rxnNotes(cytTranspIdx)));
disp('Number of cytosol transporters from the PLM:')
disp(length(rxnIdxEmptyNotes) - 2 - numAfterPLM)


%% How many intra-cellular transporters do we have?
cytIdx = find(strcmp(model.comps, 'c'));
intraTranspIdx = [];

for i = 1:length(model.rxns)
    % Get metabolite indices associated with the current reaction
    met_indices = find(model.S(:, i) ~= 0);

    % Check compartments for these metabolites
    comps = unique(model.metComps(met_indices));
    
    % If 'c' compartment is present along with other compartments
    if ~ismember(cytIdx, comps) && length(comps) > 1 && ...
        ~startsWith(model.rxns(i), 'Exch_')
        intraTranspIdx = [intraTranspIdx; i];
    end
end

disp('Number of intraorganellar transporters:')
disp(length(intraTranspIdx))

%% Functions
function comp = extractComp(metId)
    metComp = regexp(metId, '\[(.*?)\]', 'tokens');
    if ~isempty(metComp)
        comp = metComp{end}{1};
    else
        comp = '';
    end
end
