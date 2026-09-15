% Adapted from doi.org/10.1371/journal.pone.0229408 :
% - replaced the functions not published with the paper as good as possible
% - adapted to pcm
% improvements:
%   - check whether model fields are present, only then merge the fields
%   - let the user choose the symbols of chloroplast etc. compartment in
%       the cytosol model 
%   - replace chloroplast model if already present
%   - added support for merging enzyme-constrained models

function [model, missingMets] = plugAndPlay(cyt, chl, cytIsEc, bioRxnsIx)

% Merges model without chloroplast compartment with chloroplast model
% cytIsEc (optional) set to true if the cytosol model is enzyme constrained
%   (default: false)
% bioRxnsIx (optional) A logical/zero-and-one array of length of cyt.rxns
%   identifying the biomass rxns. If not given, the script will prompt the
%   user and ask whether its guess is correct.

% Additional scripts required to run plugAndPlay.m:
%   - metKEGGIDsearch.m
%   - fillKEGGIDholes.m
%   - KEGGIDfirstAid.m
%   - updateFromGrRules.m

if nargin < 3
    cytIsEc = false;
end

% request the user to check bioRxns. Thjs is important to know so we can
% decide whether a chloroplast reaction should be removed or not.
if nargin < 4
    bioRxnsIx = arrayfun(@(i) isBio(cyt, i), 1:length(cyt.rxns));
    disp(['We found ' num2str(sum(bioRxnsIx)) ' biomass reaction(s). ' ...
          'Please make sure this is correct and otherwise interrupt' ...
          newline '  and call this function with a biomassRxnIx. The ' ...
          'rxns we found:'])
    for i = 1:length(bioRxnsIx)
        if bioRxnsIx(i) ~= 0
            disp(cyt.rxns{i})
        end
    end
    question = ['Is this list of biomass reaction(s) ' ...
        'correct? (yes = 1, no = 0)'];
    bioListCorrect = input(sprintf(question));
    if ~bioListCorrect
        error(['Please call this script with an index identfying the ' ...
            'correct biomass reactions. The script excludes these from' ...
            newline '  being removed as part of the chloroplast. ' ...
            'Exchange reactions are handles separately.'])
    end
elseif ~islogical(bioRxnsIx)
    try
        bioRxnsIx = logical(bioRxnsIx);
    catch ME
        error('We were unable to convert bioRxnsIx to logical.')
    end
elseif length(bioRxnsIx) ~= length(cyt.rxns)
    error('bioRxnsIx has to have the same length as cyt.rxns')
end

%% 0 - Saving backup of models and requesting user input

fprintf('\nThis script will make changes to input models\n\n')
backup = input('Do you want to save a backup version of the models? (yes = 1, no = 0): ');
disp(' ')
if backup == 1
    cyt_filename = sprintf('%s_Backup_cytosol_model.mat', datetime("now"));
    chl_filename = sprintf('%s_Backup_chloroplast_model.mat', datetime("now"));
    save(cyt_filename, 'cyt')
    save(chl_filename, 'chl')
end
clear cyt_filename chl_filename backup

% request input for compartment symbols in the cytosol model
chlCompKeys = chl.comps;
chlModelCytSymb = 'c';
chlCompVals = chl.compNames;
cytCompSymbs = chlCompKeys;
listStr = '';
for i = 1:length(chlCompKeys)
    listStr = [listStr chlCompKeys{i} '\t' chlCompVals{i} '\n'];
end
question = ['This script will expect the following compartment symbols' ...
            ' in the cytosol model, if the compartments are present:' ...
            '\n\n' listStr '\nAre any of these compartments or symbols' ...
            ' present in the cytosol model, but with different names or'...
            ' symbols? (yes = 1, no = 0) '];
changeMetSymbols = input(sprintf(question));
if changeMetSymbols
    fprintf(['For each compartment, please input...\n' ...
        '- the compartment''s symbol of the cytosol model\n' ...
        '- any non-existing symbol if the compartment symbol is present, but not the compartment\n' ...
        '- 0 if the compartment is not present\n' ...
        '- enter if the symbol is correct.\n']);
    for i = 1:length(chlCompKeys) 
        currInput = input(sprintf(['Change ' chlCompKeys{i} '\t' ...
            chlCompVals{i} '?']), "s");
        if ~isempty(currInput)
            cytCompSymbs{i} = currInput;
        end
    end
else
    disp(['Trying to get the compartments of cytosol model from ' ...
        'cyt.comps!']);
    for i = 1:length(cytCompSymbs)
        if ~any(strcmp(cyt.comps, cytCompSymbs{i}))
            cytCompSymbs{i} = '0';
        end
    end
    if strcmp(cytCompSymbs(strcmp(chlCompKeys, chlModelCytSymb)), '0')
        error(['Could not identify the cytosol compartment in the ' ...
            'cytosol model']);
    end
end
clear currInput changeMetSymbols question listStr

% check which of the gene fields are present in the cytosol model
mergeGeneFields = true;
if ~isfield(cyt, 'grRules')
    warning(['This script will generate the fields rxnGeneMat, genes, ' ...
        'and rules from the grRules field. The cytosol model does not' ...
        ' appear to contain this field. The script will not attempt ' ...
        'to merge the gene fields']);
    mergeGeneFields = false;
end

%% 1 - Changes to chl model:
% - Deleting biomass and exchange reactions
% - Remove now-orphaned metabolites
% - Remove pseudo metabolites
% - Change compartment symbol of cytosol reactions

bioList = smatch(chl.rxnNames, 'biomass');
bioNames = chl.rxns(bioList);
for i = 1:length(bioNames)
    chl = removeRxns(chl, bioNames(i));
end
clear bioList bioNames

% Remove (exchange, import, export) reactions
exchRxnIx = find(chl.rxn2subSystem(:, strcmp(chl.subSystemNames, 'export')) | ...
    chl.rxn2subSystem(:, strcmp(chl.subSystemNames, 'import')) | ...
    chl.rxn2subSystem(:, strcmp(chl.subSystemNames, 'exchange')) | ...
    startsWith(chl.rxns, 'Ex_') | startsWith(chl.rxns, 'Im_') | ...
    startsWith(chl.rxns, 'Exch_'));
if isfield(chl, 'confidenceScores')
    chl.confidenceScores(exchRxnIx) = [];
end
toRemove = chl.rxns(exchRxnIx);
chl = removeRxns(chl, toRemove);
for rxnIx = 1:length(toRemove)
    if any(strcmp(chl.rxns, toRemove{rxnIx}))
        disp(toRemove{rxnIx})
    end
end

% remove pseudo metabolites used to constrain reactions
constMets = setdiff(chl.mets(startsWith(chl.mets, 'CONST')), ...
    'CONST_rubisco[h]');  % this ensures photorespiration
[chl, ~] = removeMetabolites(chl, constMets, true, 'exclusive');

%% 2 - Checking compartment symbol overlap
% Checking that the compartment symbols of the two models to be merged do
% not overlap

% Determining where compartment symbol is and changing format if necessary
met = char(cyt.mets(1));
cytProtIx = zeros(length(cyt.mets), 1);
if strcmp(met(end-1), '_')
    for i = 1:length(cyt.mets) % Changing format of cyt.mets vector
        currMetId = cyt.mets{i};
        if ~cytIsEc || strcmp(currMetId(end-1), '_')
            oMet = char(cyt.mets(i));
            comp = oMet(end);
            nMet = sprintf('%s[%s]', oMet(1:end-2), comp);
            cyt.mets(i) = cellstr(nMet);
        elseif cytIsEc
            % this is most likely a protein from the cyt model, give it a
            % corresponding compartment label
            cyt.mets{i} = [cyt.mets{i} '[' cytCompSymbs{strcmp(chlCompKeys, chlModelCytSymb)} ']'];
            cytProtIx(i) = 1;
        else
            disp('Warning: encountered metabolite with unknown compartment format:')
            disp(currMetId);
        end
    end
    clear oMet comp nMet
elseif ~strcmp(met(end), ']')
    fprintf('Example of metabolite ID: %s\n\n', met)
    compSymb = input('Compartment symbol in example above: ', 's');
    disp(' ')
    fullID = input('Full metabolite ID in example above: ', 's');
    disp(' ')
    pos = input('Compartment symbol situated at the beginning (1) or the end (2)? ');
    disp(' ')
    if pos == 1
        x = 1;
        while ~strcmp(met(x), compSymb)
            x = x + 1;
        end
        IDstart = strfind(met, fullID);
    elseif pos == 2
        x = length(met);
        while ~strcmp(met(x), compSymb)
            x = x - 1;
        end
        IDend = strfind(met, fullID) + length(fullID) - 1;
    end
    clear compSymb
    for i = 1:length(cyt.mets) % Changing format of cyt.mets vector
        oMet = char(cyt.mets(i));
        comp = oMet(x);
        if pos == 1
            nMet = sprintf('%s[%s]', oMet(IDstart:end), comp);
        elseif pos == 2
            nMet = sprintf('%s[%s]', oMet(1:IDend), comp);
        end
        cyt.mets(i) = cellstr(nMet);
    end
    clear oMet comp nMet IDstart IDend
end
clear met pos fullID

% Changing format of compartment symbols in cytosol model
metNameCheck = cyt.metNames{1};
if (strcmp(metNameCheck(end), ']') == 0) || ...
        (strcmp(metNameCheck(end-2), '_') > 0) || ...
        (strcmp(metNameCheck(end-1), '_') > 0)
    if strcmp(metNameCheck(end-2), '_')
        nameFlag = 3;
    elseif strcmp(metNameCheck(end-1), '_')
        nameFlag = 2;
    else
        nameFlag = 1;
    end
    for i = 1:length(cyt.metNames)
        name = char(cyt.metNames(i));
        met = char(cyt.mets(i));
        comp = regexp(met, '\[(.*?)\]', 'tokens');
        if ~isempty(comp)
            comp = comp{end}{1};
        else
            continue;
        end
        if ~endsWith(name, [' [' comp ']'])
            if nameFlag == 1
                nName = sprintf('%s [%s]', name, comp);
            elseif nameFlag == 2
                nName = sprintf('%s [%s]', name(1:end-2), comp);
            elseif nameFlag == 3
                nName = sprintf('%s [%s]', name(1:end-3), comp);
            end

            cyt.metNames(i) = cellstr(nName);
        end
    end
end
clear metNameCheck nameFlag name met comp nName

%% 2.5 remove chloroplast compartments from the cyt model

fprintf('Removing chloroplast compartments from cytosol model.\n')
rxnToRemIx = zeros(length(cyt.rxns), 1);
protectedRxns = zeros(length(cyt.rxns), 1);
protectedMets = zeros(length(cyt.mets), 1);

% neither remove cytosol, nor compartments non-existent in cyt model
%   cytCompSymbs only encompasses comps also present in pcm
compsToRemove = setdiff(unique(cytCompSymbs), ...
    [cytCompSymbs(strcmp(chlCompKeys, chlModelCytSymb)), {'0'}]);

% Ask the user whether to remove transports to a chloroplast compartment
disp(['Compartments ' strjoin(compsToRemove, ', ') ' will be removed.']);
remTrans = input(['Should transporters to non-chloroplast compartments' ...
                  ' also be removed? (yes = 1, no = 0): ']);

% gather reactions to remove
for j = 1:length(cyt.rxns)

    isBioRxn = bioRxnsIx(j);
    isExch = isExchange(cyt, j);
    currMetsIx = find(cyt.S(:, j));

    % do not remove bio reactions or import reactions (like light)
    if ~isBioRxn && ~isExch

        % iterate over metabolites in the current reaction
        allMetsInChloro = true;
        anyMetsInChloro = false;
        for k = 1:length(currMetsIx)
    
            % ignore protein metabolites
            if cytIsEc && cytProtIx(currMetsIx(k))
                continue
            end
    
            met = cyt.mets{currMetsIx(k)};
            currComp = regexp(met, '\[(.*?)\]', 'tokens');

            if ~isempty(currComp)
    
                % Is the current metabolite in a compartment to be removed?
                currCompIsChloro = ismember(currComp{end}, compsToRemove);
    
                % Record whether all or some of the metabolites of the
                % current rxn are in compartments that should be removed
                allMetsInChloro = allMetsInChloro && currCompIsChloro;
                anyMetsInChloro = anyMetsInChloro || currCompIsChloro;
            else
                allMetsInChloro = false;
            end
        end
    
        % this only involves proteins: do nothing
        if cytIsEc && allMetsInChloro && ~anyMetsInChloro
            continue
        end
    end

    protectedRxns(j) = isBioRxn || isExch || ...
        (~remTrans && anyMetsInChloro && ~allMetsInChloro);

    rxnToRemIx(j) = ~protectedRxns(j) && ...
        (allMetsInChloro || ( remTrans && anyMetsInChloro ));
end

% for reactions protected from removal, iterate a second time
% to find and rename metabolites that remain in the chloroplast
for j = 1:length(protectedRxns)
    if protectedRxns(j)
        currMetsIx = find(cyt.S(:, j));
        for k = 1:length(currMetsIx)
            metID = cyt.mets{currMetsIx(k)};
            if isInComps(metID, compsToRemove)
                % get comp
                currComp = regexp(metID, '\[(.*?)\]', 'tokens');
                currComp = currComp{end}; % cannot be empty bc isInComps==1
                currComp = currComp{1};

                % get corresponding pcm comp
                newComp = chlCompKeys{strcmp(currComp, cytCompSymbs)};

                % rename the metabolite to match pcm comp name
                cyt.mets{currMetsIx(k)} = strrep(metID, ...
                                    sprintf('[%s]', currComp), ...
                                    sprintf('[%s]', newComp));

                % add to protectedMets
                protectedMets(currMetsIx(k)) = 1;
            end
        end
    end
end

% remove comps and compNames symbols
for i = 1:length(compsToRemove)
    currCytSymb = compsToRemove{i};

    % remove from comps, compNames
    if isfield(cyt, 'compNames')
        cyt.compNames(strcmp(cyt.comps, currCytSymb)) = [];
    end
    if isfield(cyt, 'comps')
        cyt.comps(strcmp(cyt.comps, currCytSymb)) = [];
    end
end

% removeRxns makes met indices shift, therefore put IDs into metsToMatch
metsToMatch = cyt.mets(logical(protectedMets));
% and protein met IDs in protMets
protMets = cyt.mets(logical(cytProtIx));

cyt = removeRxns(cyt, cyt.rxns(logical(rxnToRemIx)));

% delete any proteins that only occur in the draw reactions now
cytProtIx = ismember(cyt.mets, protMets);
sNotZero = (cyt.S ~= 0);
cytProtsToRemoveIx = (sum(sNotZero, 2) == 1) & cytProtIx;

clear rxnToRemIx allMetsInChloro anyMetsInChloro currCompIsChloro;
clear met newComp currMetsIx currCytSymb protectedMets cytProtIx;
clear protMets sNotZero;

% for some reason there are still some orphan metabolites in the legacy
% compartments after this
cyt = removeMetabolites(cyt, cyt.mets(cytProtsToRemoveIx | ...
                                      ~any(cyt.S ~= 0, 2)), ...
                        true, 'exclusive');

% Checking that [c] now is the only overlapping compartment symbool between
% the models

% gather all compartments of cytosol model
cytComps = cell.empty(0,1); % Making empty vector for compartments of cytosol model metabolites
for i = 1:length(cyt.mets) % Filling cytosol model metabolite vector
    met = char(cyt.mets(i));
    currComp = regexp(met, '\[(.*?)\]', 'tokens');
    if strcmp(met(end), ']') && ~isempty(currComp)
        cytComps(i,1) = currComp{end}; % Take the last match
    elseif ~cytIsEc
        disp(['Warning: something must have gone wrong with ' ...
              'conversion of IDs to correct format...'])
    end
end
clear met

cytComps = unique(cytComps);

% Check if pcm compartment symbols are found amongst compartment symbols of
% cytosol model
for ChlCompIx = 1:length(chlCompKeys)
    chlCurrSymbol = chlCompKeys{ChlCompIx};
    currCompName = chlCompVals{ChlCompIx};
    if strcmp(chlCurrSymbol, chlModelCytSymb)
        continue
    end
    if any(strcmp(cytComps, chlCurrSymbol))
        rxnsInCurrComp = zeros(1, length(cyt.rxns));
        for j = 1:length(cyt.mets)
            if contains(cyt.mets{j}, ['[' chlCurrSymbol ']']) && ...
                ~any(strcmp(metsToMatch, cyt.mets{j}))
                rxnsInCurrComp = rxnsInCurrComp + (cyt.S(j, :) ~= 0);
            end
        end
        rxnsInCurrComp = rxnsInCurrComp ~= 0;
        hRxns = cyt.rxns(rxnsInCurrComp);

        % check whether the reactions are simply biomass reactions or
        % import reactions or transporters that were not removed on purpose
        harmlessRxns = 1;
        for j = 1:length(hRxns)
            currRxnIx = strcmp(cyt.rxns, char(hRxns(j)));
            if ~contains(lower(char(hRxns(j))), 'bio') && ...
                    ~isExchange(cyt, currRxnIx) && ...
                ~(~remTrans && isTransport(cyt, currRxnIx))
                harmlessRxns = 0;
            else
                % if it is harmless and being kept, we have to match up the
                % chloroplast metabolites at some point

            end
        end
        if harmlessRxns == 0
            disp(['Symbol for ' currCompName ' in chloroplast model (', ...
                chlCurrSymbol, ') is found in cytosol model']);
            disp(' ');
            newChloroplastSymb = input(['Choose new ' currCompName ' symbol for chloroplast model: '], 's');
            disp(' ')
            for i = 1:length(chl.mets) % Changing compartment symbol in chloroplast model
                met = char(chl.mets(i));
                if strcmp(met(end-1), chlCurrSymbol)
                    newMetID = sprintf('%s[%s]', met(1:end-3), newChloroplastSymb);
                    chl.mets(i) = cellstr(newMetID);
                    oldMetName = chl.metNames{i};
                    newMetName = sprintf('%s [%s]', oldMetName(1:end-4), newChloroplastSymb);
                    chl.metNames(i) = cellstr(newMetName);
                end
            end
            chlCompKeys{strcmp(chlCompKeys, chlCurrSymbol)} = newChloroplastSymb;
        end
        clear newChloroplastSymb met newMetID oldMetName newMetName currentEq hRxns harmlessRxns j rxnsInCurrComp
    end
end
clear cytComps met i ChlCompIx chlCurrSymbol currCompName remTrans;

disp(['From imports/biomass/transporters, ' ...
    num2str(length(metsToMatch)) ' metabolites need to be matched.'])


%% 3 - Making sure cytosol model contains KEGG IDs for metabolites
% convert nested cell arrays to flat cell arrays where multiple KEGG IDs
%   are separated by comma
if ~isfield(cyt, 'metKEGGID')
    cyt.metKEGGID = cell(size(cyt.rxns));
end
if iscell(cyt.metKEGGID{1})
    cyt.metKEGGID = cellfun(@ (x) strjoin(x, ','), cyt.metKEGGID, ...
        'UniformOutput', false);
end
if iscell(cyt.rxnKEGGID{1})
    cyt.rxnKEGGID = cellfun(@ (x) strjoin(x, ','), cyt.rxnKEGGID, ...
        'UniformOutput', false);
end

KEGGsPresent = ~cellfun(@isempty, cyt.metKEGGID); % Checking how many metabolites in cytosol model that has a KEGG ID
KEGGID = length(find(KEGGsPresent == 1)); % Determining amount of metabolites that do or do not have KEGG IDs
noKEGGID = length(find(KEGGsPresent == 0));

if (KEGGID / (KEGGID + noKEGGID)) < 0.6
    fprintf('KEGG IDs of metabolites are necessary to connect similar metabolites of the two models\n');
    fprintf('%i out of %i metabolites in the cytosol model are not associated with a KEGG ID.\n\n', noKEGGID, length(cyt.mets));
    fill = input('Fill in KEGG-IDs? Warning: This is a feature from the original script by Roekke et al and was not tested with the pcm. (yes (use detailed method) = 1, yes (use quick method) = 2, no = 0): ');
    disp(' ')
    if fill == 1 % Detailed fill-mode
        cyt = metKEGGIDsearch(cyt); % Attemting to find more metabolite KEGG IDs
        
        stillEmptyKEGGs = 0; % Checking if there are more KEGG ID holes to be filled
        for i = 1:length(cyt.metKEGGID)
            if isempty(cyt.metKEGGID{i})
                stillEmptyKEGGs = 1;
            end
        end
        if stillEmptyKEGGs == 1 % Filling more KEGG IDs
            fprintf('Attempting to fill in the last missing KEGG IDs. Hold on tight...\n\n')
            cyt = fillKEGGIDholes(cyt);
        end
        
        stillEmptyKEGGs = 0; % Checking if there are metabolites without KEGG ID after two rounds of filling
        for i = 1:length(cyt.metKEGGID)
            if isempty(cyt.metKEGGID{i})
                stillEmptyKEGGs = 1;
            end
        end
        
        if stillEmptyKEGGs == 0 % Celebrating :-D
            fprintf('Every metabolite now has a KEGG ID! Hooray :-D\n\n')
        end
    elseif fill == 2
        fprintf('Attempting to find KEGG-IDs only for the metabolites being transported between the chloroplast and the cytoplasm\n\n')
        cyt = KEGGIDfirstAid(cyt, chl);
    end
    if fill > 0
        backup = input('Save backup of cytosol model now that KEGG IDs have been added? (yes = 1, no = 0): ');
        disp(' ')
        if backup == 1
            filename = input('Filename? (.mat will be added automatically): ', 's');
            disp(' ')
            filename = sprintf('%s.mat', filename);
            save(filename, 'cyt');
        end
    end
end
clear i fill KEGGID KEGGsPresent noKEGGID stillEmptyKEGGs backup

%% 4 - Translating exchange metabolite of chloroplast model to cytosol namespace

fprintf('Translating exchange metabolites in chloroplast model into cytosol model namespace\n\n')

% Generate a list of cytosol metabolites in chloroplast model
cMets = find(strcmp( ...
    cellfun(@extractComp, chl.mets, 'UniformOutput', false), ...
    chlModelCytSymb));

% Identify the versions of cytosol metabolites in other compartments
compMets = zeros(length(chlCompKeys), length(cMets));
for compIx = 1:length(chlCompKeys)
    for a = 1:length(cMets)
        i = cMets(a);
        met = chl.mets{i};
        pMet = sprintf('%s[%s]', met(1:end-3), chlCompKeys{compIx});
        pMatches = smatch(chl.mets, pMet, 'exact');
        if ~isempty(pMatches) && isscalar(pMatches)
            compMets(compIx,a) = pMatches;
        end
    end
end
clear i met n pMet a compIx pMatches

% Before changing the name of chloroplast model metabolites, check whether
% any metabolites in the chloroplast of the cytosol model already overlap
% without having a common KEGG ID
for i = 1:length(metsToMatch)
    if any(strcmp(metsToMatch{i}, chl.mets))
        matchingChlMetIx = strcmp(metsToMatch{i}, chl.mets);
        cytMetIx = strcmp(metsToMatch{i}, cyt.mets);
        cytKegg = cyt.metKEGGID{cytMetIx};
        chlKegg = chl.metKEGGID{matchingChlMetIx};
        keggMatch = false;
        for singleChlKegg = strsplit(chlKegg, ',')
            keggMatch = keggMatch || contains(cytKegg, singleChlKegg);
        end
        if ~keggMatch
            disp(['Found metabolite ' metsToMatch{i} ' in both models ' ...
                'but without matching KEGG IDs.']);
            disp(['Metabolite name in chloroplast model: ' chl.metNames{matchingChlMetIx}]);
            disp(['Metabolite name in cytosol model: ' cyt.metNames{cytMetIx}]);
            disp(' Do you want to...');
            disp('(0) append _chl to the metabolite in the chloroplast model or');
            res = input('(1) match the two metabolites anyway? ');
            if res == 0  % rename in chl model
                [currMet, compPart] = strtok(chl.mets{matchingChlMetIx}, '[');
                chl.mets{matchingChlMetIx} = [currMet '_chl' compPart];

            % do nothing but add KEGG ID to make sure the chl metabolites
            % is not automatically matched to another metabolite later
            elseif res == 1  
                if ~isempty(chlKegg)
                    if ~isempty(cytKegg)
                        cyt.metKEGGID{cytMetIx} = [cytKegg ', ' chlKegg];
                    else
                        cyt.metKEGGID{cytMetIx} = chlKegg;
                    end
                end
            end
        end
    end
end

% Make a cell array for chl model c-mets not present in cyt model
missingMets = cell.empty(0,2); 
n = 1;
cytCytosolSymb = cytCompSymbs{strcmp(chlCompKeys, chlModelCytSymb)};

% Find which KEGGIDs are duplicated in each model. These cannot be matched
% automatically
chlDupKEGG = getDupKEGGs(chl);

% Search for cMets in cytosol model, translating cMets and hMets into
% cytosol namespace
changedMetIDs = {};
cytComps = cellfun(@extractComp, cyt.mets, 'UniformOutput', false);
a = 1; % Solely for numbering metabolites in the printed output
chlMetsProcessed = zeros(size(chl.mets));
for cMetsIx = 1:length(cMets)
    i = cMets(cMetsIx);
    if chlMetsProcessed(i) == 1
        continue;
    end

    metName = chl.metNames{i};
    KEGGID = strtrim(chl.metKEGGID{i});
    chlMet = chl.mets{i};
    chlMetComp = extractComp(chlMet);

    % Find matches in cyt based on KEGG ID and compartment
    metMatch = findMatches(cyt, strsplit(KEGGID, ','), cytComps, ...
        cytCompSymbs(strcmp(chlCompKeys, chlMetComp)));

    if isempty(metMatch)
        % Met might not be present in cytosol model
        [chl, missingMets, n] = processMissingMet(chl, missingMets, ...
                    i, n, a, length(cMets), cytCytosolSymb);
    else

        if ~chlDupKEGG(i)
            % Normal metabolite without duplications

            if isscalar(metMatch)
                cytCmet = metMatch;
            else
                fprintf(['* Altering namespace for chloroplast ' ...
                    'metabolite %s\n\n  Metabolite hits in cytosol' ...
                    ' model:\n'], chlMet)
                disp([transpose(num2cell(1:length(metMatch))), ...
                    cyt.mets(metMatch)])
                correctMet = input(['Index of correct metabolite ' ...
                    '(if correct metabolite is not present, press 0): ']);
                disp(' ')
                if correctMet > 0
                    cytCmet = metMatch(correctMet);
                elseif correctMet == 0
                    cytCmet = 0;
                end
            end
    
            % Changing name of cytosol metabolite and corresponding plastid
            % metabolite(s) in chloroplast model
            if cytCmet == 0
                % Metabolite is not present in cytosol compartment of
                % cytosol model
                [chl, missingMets, n] = processMissingMet(chl, ...
                    missingMets, i, n, a, length(cMets), cytCytosolSymb);
            else
                [chl, changedMetIDs] = processMetMatch(chl, cyt, i, ...
                    cytCmet, a, cMetsIx, length(cMets), cytCytosolSymb, ...
                    chlCompKeys, chlModelCytSymb, compMets, ...
                    changedMetIDs, cytCompSymbs);
            end
            chlMetsProcessed(i) = 1;
            a = a + 1;

        else
            % if the current chl metabolite has duplicated KEGGID, ask for
            % user input

            % Search for all chl metabolites of the current set of KEGG IDs
            currChlMets = zeros(size(chl.mets)); % will be find() later
            currChlMets(i) = 1;

            % Start the search with the current chl KEGG IDs
            searchSet = strsplit(KEGGID, ',');
            finishedSet = {};
        
            while ~isempty(searchSet)
                % Pop the first KEGG ID from the search set for processing
                currentKEGGID = searchSet{1};
                searchSet(1) = [];
                finishedSet = [finishedSet, currentKEGGID];
        
                % Find chl metabolites matching the KEGG ID and compartment
                for j = 1:length(chl.metKEGGID)
                    chlKeggIds = strsplit(strtrim(chl.metKEGGID{j}), ',');
                    if ismember(currentKEGGID, chlKeggIds) && ...
                            strcmp(extractComp(chl.mets{j}), chlMetComp)
                        currChlMets(j) = 1;
        
                        % Add new KEGG IDs to the search set if not already finished
                        for keggId = chlKeggIds
                            if ~isempty(keggId) && ...
                                ~ismember(keggId, finishedSet) && ...
                                ~ismember(keggId, searchSet)
                                searchSet = [searchSet, keggId];
                            end
                        end
                    end
                end
            end
            currChlMets = find(currChlMets);
            
            % Find matches in cyt
            currCytMets = findMatches(cyt, finishedSet, cytComps, ...
                cytCompSymbs(strcmp(chlCompKeys, chlMetComp)));

            if ~isempty(currCytMets)
                % Communicate with the user
                fprintf('* Metabolite %i of %i:\n', a, length(cMets))
                fprintf(['  Encountered metabolite %s with non-unique ' ...
                        'KEGG ID %s in the chloroplast model.\n'], ...
                        metName, KEGGID)
                fprintf('  Metabolites in chl and cyt models:\n')
                
                % Pad shorter columns with empty cells
                maxLength = max(length(currChlMets), length(currCytMets));
                chlPad = repmat({''}, maxLength - length(currChlMets), 1);
                cytPad = repmat({''}, maxLength - length(currCytMets), 1);
                disp([[transpose(num2cell(1:length(currChlMets))); chlPad], ...
                      [chl.mets(currChlMets); chlPad], ...
                      [transpose(num2cell(1:length(currCytMets))); cytPad], ...
                      [cyt.mets(currCytMets); cytPad]]);

                fprintf(['  For each metabolite in the chl model, please ' ...
                    'type the index of the corresponding\n    metabolite ' ...
                    'in the cyt model (if correct metabolite is not ' ...
                    'present, press 0):\n'])
            end

            % for each chl metabolite:
            % ask the user to pick the corresponding cyt metabolite or 0
            for chlIx = 1:length(currChlMets)
                if ~isempty(currCytMets)
                    correctMet = input(['Index of correct metabolite in ' ...
                        'cyt model for (' num2str(chlIx) ') ' ...
                        chl.mets{currChlMets(chlIx)} ': ']);
                    if correctMet > 0
                        cytCmet = metMatch(correctMet);
                    elseif correctMet == 0
                        cytCmet = 0;
                    end
                else
                    cytCmet = 0;
                end

                % Change names and IDs of cytosol and corresponding plastid
                % metabolites in chl model
                if cytCmet == 0
                    [chl, missingMets, n] = processMissingMet(chl, ...
                        missingMets, currChlMets(chlIx), n, a, ...
                        length(cMets), cytCytosolSymb);
                else
                    [chl, changedMetIDs] = processMetMatch(chl, cyt, ...
                        currChlMets(chlIx), cytCmet, a, cMetsIx, ...
                        length(cMets), cytCytosolSymb, chlCompKeys, ...
                        chlModelCytSymb, compMets, changedMetIDs, ...
                        cytCompSymbs);
                end
                chlMetsProcessed(currChlMets(chlIx)) = 1;
                a = a + 1;
            end
        end
    end
end
clear cytCmet cytMet a i j m n KEGGID chlMet metMatch metName constName
clear refMet refRoot checkMet checkRoot correctMet compMets cMets cMetsIx
clear chlMetsProcessed chlDupKEGG

% After changing the name of chl model metabolites, check whether any mets
% in metToMatch are left unmatched
disp('Checking which metabolites are not ready to be mapped automatically')
chlComps = cellfun(@extractComp, chl.mets, 'UniformOutput', false);
cytDupKEGG = getDupKEGGs(cyt);
for i = 1:length(metsToMatch)
    if ~any(strcmp(metsToMatch{i}, chl.mets))
        cytMetIx = strcmp(metsToMatch{i}, cyt.mets);
        cytKegg = cyt.metKEGGID{cytMetIx};

        % match KEGG IDs provided compartments match up
        matchingChlMets = findMatches(chl, strsplit(cytKegg, ','), ...
            chlComps, extractComp(metsToMatch{i}));

        chlMetIx = 0;
        % if we find one match, we have a match unless the chl met was 
        % changed before or the cyt met is one of many with that KEGG ID
        if isscalar(matchingChlMets)
            chlMetIx = matchingChlMets(1);
            if any(strcmp(changedMetIDs, chl.mets{chlMetIx}))
                % was matched before; put cytMet into missing Mets
                missingMets(end + 1, 1) = cellstr(cyt.metNames{cytMetIx});
                missingMets(end, 2) = cellstr([cyt.mets{cytMetIx} '; ' cytKegg]);
                fprintf('* Chloroplast metabolite %i of %i:\n', i, length(metsToMatch));
                fprintf('  Metabolite %s was kept in cytosol model, but its match %s was changed due to a match in the cytosol before\n', ...
                    cyt.metNames{cytMetIx}, chl.metNames{chlMetIx});
                fprintf('  Metabolite has been added to list of missing metabolites\n\n')
                chlMetIx = -1;
            elseif cytDupKEGG(cytMetIx)
                % multiple cyt mets in that compartment have that KEGG ID
                currMatchingCytMets = findMatches(cyt, strsplit(cytKegg, ','), ...
                    cytComps, extractComp(metsToMatch{i}));

                % ask user with all possible cyt mets
                fprintf('* Chloroplast metabolite %i of %i:\n', i, length(metsToMatch));
                fprintf('  Trying to match following chloroplast metabolites from cyt model to chl model.\n\n')
                disp([transpose(num2cell(1:length(currMatchingCytMets))), cyt.mets(currMatchingCytMets)])
                fprintf('  For their KEGG ID %s, found %s as hit.\n', cytKegg, chl.mets{matchingChlMets})
                correctMet = input('Index of correct cyt metabolite (if correct metabolite is not present, press 0): ');
                disp(' ')

                % Match the chosen metabolite
                if correctMet > 0
                   oldPid = chl.mets{chlMetIx};
                   chl.mets{chlMetIx} = cyt.mets{currMatchingCytMets(correctMet)};
                   fprintf('  Chloroplast metabolite %i of %i: Metabolite ID changed from %s to %s\n', ...
                            i, length(metsToMatch), oldPid, chl.mets{chlMetIx})
                end

                % If the chosen metabolite is not the one of the current 
                % iteration, add that one to missing metabolites
                if correctMet == 0 || ...
                        currMatchingCytMets(correctMet) ~= find(cytMetIx)
                    missingMets(end + 1, 1) = cellstr(cyt.metNames{cytMetIx});
                    missingMets(end, 2) = cellstr([cyt.mets{cytMetIx} '; ' cytKegg]);
                    fprintf('  Metabolite %s has been added to list of missing metabolites\n\n', cyt.mets{cytMetIx})
                end
                chlMetIx = -1;  % Do not do anything more
            end

        % if we find multiple matches, let user decide
        elseif length(matchingChlMets) > 1
            fprintf('* Trying to match chloroplast metabolite %s in cyt model to chl model.\n\n', cyt.mets{cytMetIx})
            fprintf('  Metabolite name: %s\n', cyt.metNames{cytMetIx})
            fprintf('  Metabolite hits in chl model for KEGG ID %s:\n', cytKegg)
            disp([transpose(num2cell(1:length(matchingChlMets))), chl.mets(matchingChlMets), chl.metNames(matchingChlMets)])
            correctMet = input('Index of correct metabolite (if correct metabolite is not present, press 0): ');
            disp(' ')
            if correctMet > 0
                chlMetIx = matchingChlMets(correctMet);
            elseif correctMet == 0
                chlMetIx = 0;
            else
                disp('Did not understand input, interpret as 0')
                chlMetIx = 0;
            end
        end

        if chlMetIx > 0
           oldPid = chl.mets{chlMetIx};
           chl.mets{chlMetIx} = metsToMatch{i};
           fprintf('  Chloroplast metabolite %i of %i: Metabolite ID changed from %s to %s\n', ...
                    i, length(metsToMatch), oldPid, char(chl.mets(chlMetIx)))

        % if we find no matches, or something else went wrong, put into
        % missingMets
        elseif chlMetIx == 0
            missingMets(end + 1, 1) = cellstr(cyt.metNames{cytMetIx});
            missingMets(end, 2) = cellstr([cyt.mets{cytMetIx} '; ' cytKegg]);
            fprintf('* Chloroplast metabolite %i of %i:\n', i, length(metsToMatch));
            fprintf('  Metabolite %s was kept in cytosol model, but does not match to any metabolite in the chl model\n', cyt.metNames{cytMetIx});
            fprintf('  Metabolite has been added to list of missing metabolites\n\n')
        end
    end
end
clear i chlMetIx cytMetIx cytKegg

%% 5 - Merging models

fprintf('Merging models. Please wait...\n\n')
model = cyt;
rxnNamesModel = ReactionNames(model);
rxnNamesChl = ReactionNames(chl);

% Check if SLIMES are present and prompt user for behavior
numSLIMEs = sum(startsWith(cyt.rxns, 'SLIMEr_'));
slimeOption = 3;
if numSLIMEs > 10
    fprintf(['Your model also uses SLIME reactions (n = ' ...
        num2str(numSLIMEs) ') that might clash with the SLIME ' ...
        'reactions\nof the chloroplast model (n = ' ...
        num2str(sum(startsWith(chl.rxns, 'SLIMEr_'))) ...
        '). Please choose an option how to proceed:\n']);
    fprintf('(0) Let cyt SLIME overwrite chloroplast SLIME reactions\n');
    fprintf('(1) Let chloroplast SLIME overwrite cyt SLIME reactions\n');
    fprintf(['(2) Always keep both and append _chl to the chloroplast ' ...
        'SLIME reaction IDs (Recommended)\n']);
    fprintf('(3) Decide for each reaction separately (Default)\n');

    slimeOption = input('Option for handling SLIME reactions (0|1|2|3): ');
    if slimeOption > 3 || slimeOption < 0
        fprintf('Did not understand input, using default option (3).');
        slimeOption = 3;
    end

    fprintf(' ')
end

% if true, cyt uses the subsystems format with cells, otherwise the one
% with chars
cytSubSystemsFormat = any(cellfun(@iscell, cyt.subSystems));
for j = 1:length(chl.rxns) % Looping through chloroplast model reactions
    n = length(model.rxns) + 1; % Index of next reaction from chl model to be added to model
    % change name of prot_pool_exchange of chloroplast automatically 
    % (S is changed below)
    if cytIsEc && strcmp(chl.rxns{j}, 'prot_pool_exchange')
        chl.rxns{j} = 'prot_pool_chl';
    end
    % when transporters are kept, there might be a conflict of identical
    % rxns but with different stoichiometry
    if any(strcmp(model.rxns, chl.rxns(j)))
        modelRxnIx = strcmp(model.rxns, chl.rxns(j));

        % if stoichiometries are the same, just keep original
        if areRxnsEquivalent(rxnNamesModel{modelRxnIx}, rxnNamesChl{j})
            fprintf(['Rxn %s already exists with identical ' ...
                'stoichiometry. Retaining original reaction\n'], ...
                chl.rxns{j});
            continue;
        end

        % we would have to ask the user but maybe they set a global option
        if startsWith(chl.rxns{j}, 'SLIMEr_') && (slimeOption ~= 3)
            res = slimeOption;
        else
            % otherwise ask the user
            fprintf(['Rxn %s (%i / %i) already exists. Please choose which ' ...
                'version to keep (0 is from cyt, 1 from chl model):\n'], chl.rxns{j}, j, length(chl.rxns));
            disp(['(0) ' rxnNamesModel{modelRxnIx}]);
            disp(['(1) ' rxnNamesChl{j}]);
            disp('(2) Keep both, appending _chl to chl rxn ID');
            res = input('Index of reaction to keep (0|1|2): ');
            disp(' ')
        end

        if res == 0
            continue
        elseif res == 1
            model = removeRxns(model, model.rxns(modelRxnIx));
            rxnNamesModel(modelRxnIx) = [];
            n = length(model.rxns) + 1;
        elseif res == 2
            chl.rxns{j} = [chl.rxns{j} '_chl'];
        else
            disp('Input not recognized; keeping the original rxn (0)');
        end
    end

    % rxns
    model.rxns(n) = chl.rxns(j);
    % rxnName
    model.rxnNames(n) = chl.rxnNames(j);
    % subSystems
    if cytSubSystemsFormat
        model.subSystems{n} = chl.subSystemNames(logical(chl.rxn2subSystem(j, :)));
    else
        model.subSystems(n) = chl.subSystems(j);
    end
    % ub
    model.ub(n) = chl.ub(j);
    % lb
    model.lb(n) = chl.lb(j);
    % c
    model.c(n) = 0;
    % c
    if isfield(model, 'C')
        model.C(:, n) = zeros(size(model.C, 1), 1);
    end
    % grRules
    model.grRules(n) = chl.grRules(j);
    % EC numbers
    if isfield(model, 'rxnECNumbers') && isfield(chl, 'rxnEC')
        model.rxnECNumbers(n) = chl.rxnEC(j);
    end


    % handle rxn fields present in both chl and cyt, or only in cyt
    rxnFieldNames = {'rev', 'confidenceScores', 'rxnConfidenceScores', ...
                     'rules', 'rxnReferences', 'rxnNotes', 'rxnKEGGID', ...
                     'rxnisECID', 'rxnisbiggID', 'rxnismetacycID', ...
                     'rxnisseedID'};
    for i = 1:length(rxnFieldNames)
        rxnFieldName = rxnFieldNames{i};
        if isfield(model, rxnFieldName) && isfield(chl, rxnFieldName)
            model.(rxnFieldName)(n) = chl.(rxnFieldName)(j);
        elseif isfield(model, rxnFieldName)
            % figure out type of array/cell content and fill with empty
            if iscell(model.(rxnFieldName)(end)) && ...
                    ischar(model.(rxnFieldName){end})
                model.(rxnFieldName)(n) = {char.empty};
            elseif isnumeric(model.(rxnFieldName)(end))
                model.(rxnFieldName)(n) = nan;
            else
                warning('Problem with reaction fields while merging:');
                disp(['Reaction field ' rxnFieldName ' is not present' ...
                    ' in chl and its type could not be inferred to pad' ...
                    ' with empty'])
            end
        end
    end
    
    % Handling metabolites involved in reaction
    % met index in chl, S-coeff in chl, 
    % met index in model (for now set to zero),
    % metabolite found in model already?
    metMatrix = [find(chl.S(:,j) ~= 0), nonzeros(chl.S(:,j)), ...
        zeros(length(nonzeros(chl.S(:,j))),1), ...
        zeros(length(nonzeros(chl.S(:,j))),1)];
    metIDs = chl.mets(metMatrix(:,1));
    m = length(model.mets) + 1;
    for i = 1:length(metIDs)
        met = char(metIDs(i));
        modelMetIndex = smatch(model.mets, met, 'exact');
        if ~isempty(modelMetIndex)
            if isscalar(modelMetIndex)
                metMatrix(i,3) = modelMetIndex;
                metMatrix(i,4) = 1;
            else
                disp(['Warning: Metabolite %s is found more than once' ...
                    ' in metabolite vector of model'], met)
                disp(' ')
            end
        elseif isempty(modelMetIndex)
            metMatrix(i,3) = m;
            m = m + 1;
        end
    end
    
    metFieldNames = {'metFormulas', 'metChEBIID', ...
                'metPubChemID', 'metInChIString', 'metCharge', ...
                'metCharges', 'metKEGGID', 'metNotes', 'metSEEDID', ...
                'metisinchikeyID', 'metLIPIDMAPSID', 'metisInChIKeyID', ...
                'metisSMILESID', 'metisbiggID', 'metismetacycID', ...
                'metismetanetxID', 'metispubchemID', 'metisseedID'};
    for i = 1:size(metMatrix,1)
        h = metMatrix(i,1); % Metabolite index in chl model
        c = metMatrix(i,3); % Metabolite index in model
        if metMatrix(i,4) == 0  % metabolite is new in the model
            % mets
            model.mets(c) = chl.mets(h);
            % metNames
            model.metNames(c) = chl.metNames(h);
            % all other fields
            for metFieldIx = 1:length(metFieldNames)
                fieldName = metFieldNames{metFieldIx};
                if isfield(model, fieldName)
                    if isfield(chl, fieldName)
                        model.(fieldName)(c) = chl.(fieldName)(h);
                    else
                        model.(fieldName)(c) = {''};
                    end
                end
            end
        elseif metMatrix(i,4) == 1 % metabolite exists in the model

            for metFieldIx = 1:length(metFieldNames)
                fieldName = metFieldNames{metFieldIx};

                % if field exists in target model and the current entry is
                % empty, replace with chl entry. Otherwise do nothing
                if isfield(model, fieldName) && ...
        (iscell(model.(fieldName)) && isempty(model.(fieldName){c}) || ...
        (~iscell(model.(fieldName)) && isempty(model.(fieldName)(c))))
                    if isfield(chl, fieldName)
                        model.(fieldName)(c) = chl.(fieldName)(h);
                    else
                        model.(fieldName)(c) = {''};
                    end
                end
            end
        end
        % S
        model.S(c,n) = metMatrix(i,2);
    end
end

% handle chloroplast protein pool
chlPoolIdx = find(strcmp(model.rxns, 'prot_pool_chl'));
if cytIsEc && ~isempty(chlPoolIdx)
    protPoolMetIx = find(strcmp(model.mets, 'prot_pool[c]'));
    if isempty(protPoolMetIx)
        disp('Warning: Did not find protein pool metabolite of cyt model (expected prot_pool[c])')
        disp('         Please add the cytosol protein pool metabolite to rxn prot_pool_chl with a factor of 1')
    else
        model.S(protPoolMetIx, chlPoolIdx) = 1;
    end
end

% Handle genes
if mergeGeneFields
    model = updateFromGrRules(model);
end

% Add new description
if isfield(cyt, 'description')
    cytDescr = cyt.description;
elseif isfield(cyt, 'name')
    cytDescr = cyt.name;
else
    cytDescr = 'cytosolModel';
end
model.description = sprintf('Merged model made from %s and %s', ...
    cytDescr, chl.description);

clear i j n metMatrix metIDs m met h c rxGenes rxGeneNames x b yPos
clear currentSize cytDescription chlDescription borderRxns modelMetIndex

% update comps, compNames
if isfield(model, 'comps') && isfield(model, 'compNames')

    for i = 1:length(model.mets)
        met = model.mets{i};
        currComp = regexp(met, '\[(.*?)\]', 'tokens');
        if ~isempty(currComp)
            currComp = currComp{end}{1};
            if ~any(strcmp(currComp, model.comps))
                if ~any(strcmp(chlCompKeys, currComp))
                    warning(['Compartment ' currComp ' found in mets' ...
                        ' but not in model.comps nor chlCompKeys']);
                else
                    model.compNames(end + 1) = ...
                        chlCompVals(strcmp(chlCompKeys, currComp));
                    model.comps{end + 1} = currComp;
                end
            end
        end
    end
end

% handle LP fields with incorrect sizes
if isfield(model, 'b')
model.b = zeros(length(model.mets), 1);
end
if isfield(model, 'csense')
model.csense = repmat('E', length(model.mets), 1);
end
if isfield(model, 'C')
model.C = [model.C, zeros(size(model.C, 1), ...
    length(model.rxns) - size(model.C, 2))];
end

fprintf(['Congratulations, you successfully merged the cytosol ' ...
    'model and the chloroplast model!\n\nPlease run verifyModel(model)' ...
    ' to check if this script missed any field. Inspect the ' ...
    'missingMets variable and possibly call mergeMets to merge any ' ...
    'of the missing metabolites.\n\n'])
end

%% Functions not published with original Roekke et al paper
function out = smatch(str, pat, option)
    % Default to 'contains' functionality
    if nargin < 3 || ~strcmp(option, 'exact')
        out = find(contains(str, pat));
    else
        % Implement 'exact' functionality using strcmp
        out = find(strcmp(str, pat));
    end
end

function out = ReactionNames(model)
    % input: a model
    % output: cell array, one entry per reaction in model, string is like:
    % <reaction ID> : <equation of the reaction>
    numReactions = length(model.rxns);
    out = cell(numReactions, 1);
    modelWithEquations = buildRxnEquations(model);
    
    for i = 1:numReactions
        reactionID = model.rxns{i};
        equation = modelWithEquations.rxnEquations{i};
        out{i} = sprintf('%s : %s', reactionID, equation);
    end
end

%% additional functions
function isExch = isExchange(model, rxnIx)
    % checks whether the given reaction is an exchange, by checking
    % whether all metabolites are consumed, or all metabolites are
    % produced

    currFactors = model.S(:, rxnIx);
    isExch = all(currFactors >= 0) || all(currFactors <= 0);
end

function isBioRxn = isBio(model, rxnIx)
    % checks whether the given reaction is a biomass reaction. This is very
    % coarse an the user is asked to correct this in the beginning of the
    % script
    isBioRxn = contains(model.rxns{rxnIx}, 'bio', 'IgnoreCase', true) || ...
        contains(model.rxnNames{rxnIx}, 'bio', 'IgnoreCase', true);
end

function isTransp = isTransport(model, rxnIx)
    % checks whether the given reaction is an exchange, by checking
    % whether all metabolites are consumed, or all metabolites are
    % produced

    currMetsIx = find(model.S(:, rxnIx));
    currMets = model.mets(currMetsIx);
    comps = {};
    for i = 1:length(currMets)
        currComp = regexp(currMets{i}, '\[(.*?)\]', 'tokens');
        currComp = currComp{end};
        comps = [comps; currComp];
    end
    if length(unique(comps)) > 1
        isTransp = true;
    elseif isscalar(unique(comps))
        isTransp = false;
    end
end

function inComps = isInComps(metId, comps)
    metComp = extractComp(metId);
    inComps = false;
    if ~isempty(metComp)
        for compIx = 1:length(comps)
            if strcmp(metComp, comps{compIx})
                inComps = true;
                break
            end
        end
    end
end

function comp = extractComp(metId)
    metComp = regexp(metId, '\[(.*?)\]', 'tokens');
    if ~isempty(metComp)
        comp = metComp{end}{1};
    else
        comp = '';
    end
end

function res = areRxnsEquivalent(rxnEq1, rxnEq2)
    % Check if two reaction strings are equivalent, irrespective of
    % metabolite order

    % first check rxn names
    nameEq1 = split(rxnEq1, ':');
    nameEq2 = split(rxnEq2, ':');
    if ~strcmp(nameEq1{1}, nameEq2{1})
        res = false;
        return
    end

    % then check equations
    rxnEq1 = strtrim(nameEq1{2});
    rxnEq2 = strtrim(nameEq2{2});

    % split into substrates and products, if reversibility is the same
    if contains(rxnEq1, '<==>')
        parts1 = split(rxnEq1, '<==>');
        if contains(rxnEq2, '<==>')
            parts2 = split(rxnEq2, '<==>');
        else
            res = false;
            return
        end
    else
        parts1 = split(rxnEq1, '-->');
        if contains(rxnEq2, '<==>')
            res = false;
            return
        else
            parts2 = split(rxnEq2, '-->');
        end
    end
   
    % one part for substrates and one part for products
    if numel(parts1) ~= 2 || numel(parts2) ~= 2
        error('Invalid reaction format');
    end

    % Sort substrates and products for comparison
    substrates1 = sortMetabolites(parts1{1});
    products1 = sortMetabolites(parts1{2});
    substrates2 = sortMetabolites(parts2{1});
    products2 = sortMetabolites(parts2{2});

    % Compare sorted substrates and products
    res = isequal(substrates1, substrates2) && ...
            isequal(products1, products2);
end

function metabolites = sortMetabolites(part)
    % Split and sort metabolites within each part
    metabolites = sort(strtrim(split(strtrim(part), '+')));
end

function dupIx = getDupKEGGs(model)
    % Find the metabolites of the model that have non-unique KEGGIDs in
    % the model among metabolites of the same compartment
    metKEGGIDs = model.metKEGGID;
    mets = model.mets;
    
    dupIx = zeros(size(metKEGGIDs));
    keggMap = containers.Map('KeyType', 'char', 'ValueType', 'any');
    
    for i = 1:length(metKEGGIDs)
        currIDs = strsplit(metKEGGIDs{i}, ',');
        compI = extractComp(mets{i});
        if isempty(compI)
            continue
        end
        
        % Iterate over each individual ID
        for j = 1:length(currIDs)
            currID = strtrim(currIDs{j});
            if isempty(currID)
                continue;
            end
            
            % Check if the ID is already in the map
            if isKey(keggMap, currID)
                prevIx = keggMap(currID);
                
                % Check if any of the existing indices have matching comp
                for k = prevIx
                    compK = extractComp(mets{k});
                    if ~isempty(compK) && strcmp(compI, compK)
                        % If they match, flag both indices
                        dupIx(k) = 1;
                        dupIx(i) = 1;
                        break;
                    end
                end
    
                % Add to the map
                keggMap(currID) = [keggMap(currID), i];
            else
                keggMap(currID) = i;
            end
        end
    end

    dupIx = logical(dupIx);
end

% Functions for metabolite matching
% cytCompSymbs(strcmp(chlCompKeys, chlMetComp))
function currCytMets = findMatches(targetModel, keggIdSet, targetComps, ...
    searchedCompTargetNamespace)
    % Search for corresponding metabolites in targetModel and filter by 
    % compartment

    currCytMets = zeros(size(targetModel.mets)); % will be find() later
    for keggIdIx = 1:length(keggIdSet)
        keggId = keggIdSet{keggIdIx};
        if ~isempty(keggId)
            currCytMets(smatch(targetModel.metKEGGID, keggId)) = 1;
        end
    end

    % Filter for compartment
    currCytMets = find(currCytMets & strcmp(targetComps, ...
        searchedCompTargetNamespace));
end

function [chl, missingMets, n] = processMissingMet(chl, missingMets, i, n, a, numMets, cytCytosolSymb)

    metName = chl.metNames{i};
    KEGGID = strtrim(char(chl.metKEGGID(i)));
    chlMet = chl.mets{i};

    missingMets(n,1) = cellstr(metName);
    missingMets(n,2) = cellstr(KEGGID);
    n = n + 1;
    fprintf('* Metabolite %i of %i:\n', a, numMets)
    fprintf('  Metabolite %s (%s) was marked as missing in cytosol model\n', metName(1:end-3), chlMet)
    fprintf('  Metabolite has been added to list of missing metabolites\n\n')

    % Still change the compartment identifier
    chl.mets{i} = sprintf('%s[%s]', chlMet(1:end-3), cytCytosolSymb);
end

function [chl, changedMetIDs] = processMetMatch(chl, cyt, i, cytCmet, a, cMetsIx, numMets, ...
    cytCytosolSymb, chlCompKeys, chlModelCytSymb, compMets, changedMetIDs, cytCompSymbs)
    chlMet = chl.mets{i};

    % Changing cytosol metabolite ID in chl model
    metIDchange = false;
    oldCid = char(chl.mets(i));
    if ~strcmp(cyt.mets(cytCmet), chl.mets(i)) && ...
            any(strcmp(chl.mets, char(cyt.mets(cytCmet))))
        fprintf('* Metabolite %i of %i: Cannot change ID of metabolite %s to %s. %s found in chl model already. Names will still be changed.\n', ...
            a, numMets, oldCid, cyt.mets{cytCmet}, cyt.mets{cytCmet})

        % Still change the compartment identifier
        chl.mets{i} = sprintf('%s[%s]', chlMet(1:end-3), cytCytosolSymb);
    else
        metIDchange = true;
        chl.mets(i) = cyt.mets(cytCmet);
        newCid = chl.mets{i};
        fprintf('* Metabolite %i of %i: Metabolite ID changed from %s to %s\n', ...
            a, numMets, oldCid, newCid)
        changedMetIDs = [changedMetIDs, newCid];
    end

    % Changing cytosol metabolite name
    oldCname = chl.metNames{i};
    chl.metNames(i) = cyt.metNames(cytCmet);
    fprintf('  Metabolite %i of %i: Metabolite name changed from %s to %s\n', ...
        a, numMets, oldCname, chl.metNames{i})
    cytCytosol = ['[' cytCompSymbs{strcmp(chlCompKeys, chlModelCytSymb)} ']'];
    newCname = chl.metNames{i};

    % Changing other compartments metabolite ID (if appropriate) and name
    for compSymbolIx = 1:length(chlCompKeys)
        currCompSymbol = chlCompKeys{compSymbolIx};
        if strcmp(currCompSymbol, chlModelCytSymb) % cytosol was done above
            continue
        end
        currCompMetIx = compMets(strcmp(chlCompKeys, ...
                                 currCompSymbol), cMetsIx);
        if currCompMetIx == 0
            continue
        end
        
        % Change metabolite IDs of metabolites in other compartments
        if metIDchange
            oldPid = chl.mets{currCompMetIx};
            newPid = sprintf('%s[%s]', strtok(newCid, '['),currCompSymbol);
            chl.mets(currCompMetIx) = cellstr(newPid);
            fprintf('  Metabolite %i of %i: Metabolite ID changed from %s to %s\n', ...
                a, numMets, oldPid, chl.mets{currCompMetIx})
            changedMetIDs = [changedMetIDs, newPid];
        end

        % Change metabolite names of metabolites in other compartments
        oldPname = chl.metNames{currCompMetIx};
        % From newCname, remove [c] (or [cy] or whatever)
        if endsWith(strip(newCname), cytCytosol)
            newCname = strip(newCname);
            newCname = strip(newCname(1:end-length(cytCytosol)));
        end
        newPname = sprintf('%s [%s]', newCname, currCompSymbol);
        chl.metNames(currCompMetIx) = cellstr(newPname);
        fprintf('  Metabolite %i of %i: Metabolite name changed from %s to %s\n', ...
            a, numMets, oldPname, chl.metNames{currCompMetIx})
    end
    disp(' ')
end
