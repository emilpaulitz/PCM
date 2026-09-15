%% Add KEGG reaction to model
function [model, addedRxnId, errorMsg] = addKEGGReactionNew(model, ...
    KEGG_ID, gpr, rxnNote, subSystems)
    if any(contains(model.rxns,KEGG_ID))
        disp(['Reaction ' KEGG_ID ' already in the model'])
        addedRxnId = char.empty;
        errorMsg = '';
        return 
    end

    % do request
    [r_entry, success] = ...
        tryRequest('https://rest.kegg.jp/get/reaction:', KEGG_ID);
    if ~success
        addedRxnId = char.empty;
        errorMsg = 'Failed initial webrequest';
        return
    end

    reaction=convertStringsToChars(regexp(string(r_entry), ...
        "(?<=EQUATION\s+)\S.+?(?=\n)",'match','once'));

    % catch glycan reactions
    if ~isempty(regexp(reaction, '[G][0-9]+', 'match', 'once'))
        newId = convertStringsToChars(regexp(string(r_entry), ...
            "(?<=REMARK\s+Same as:\s+)\S.+?(?=\n)",'match','once'));
        if any(contains(model.rxns, newId))
            disp([KEGG_ID ' is glycan and its alternative is in model:']);
            disp(model.rxns(contains(model.rxns, newId)));
            addedRxnId = char.empty;
            errorMsg = char.empty;
            return
        else
            disp([KEGG_ID ' is glycan, adding its alternative: ' newId]);
            KEGG_ID = newId;

            [r_entry, success] = ...
                tryRequest('https://rest.kegg.jp/get/reaction:', KEGG_ID);
            if ~success
                addedRxnId = char.empty;
                errorMsg = 'Failed webrequest for alternative';
                return
            end
            reaction=convertStringsToChars(regexp(string(r_entry), ...
                    "(?<=EQUATION\s+)\S.+?(?=\n)",'match','once'));
        end
    end

    % catch reactions with the same compound on both sides
    mets = regexp(reaction, '[C][0-9]+', 'match');
    if length(unique(mets)) ~= length(mets)
        disp(['Error: the same compound occurs multiple times ' KEGG_ID])
        addedRxnId = char.empty;
        errorMsg = 'The same compound occurs multiple times';
        return
    end

    % read other properties from the result
    r_name = convertStringsToChars(regexp(string(r_entry), ...
        "(?<=NAME\s+)\S.+?(?=\n)",'match','once'));
    ecs = strsplit(convertStringsToChars(regexp(string(r_entry), ...
        "(?<=ENZYME\s+)\S.+?(?=\n)",'match','once')));
    ecs = remGeneralECs(ecs);
    
    % put mets into compartments
    [rxnEquation, newMets, chosenComp, errorMsg] = modifyForModel(reaction, model);

    % request information on new metabolites
    for i = 1: length(newMets)
        met = newMets{i};
        res = readAndParseMet(met);
        if isempty(fieldnames(res))
            disp(['Error: no fields found for met ' met])
            addedRxnId = char.empty;
            errorMsg = [errorMsg 'Failed fetching metabolite details for ' met '; '];
            return
        end
        if contains(res.formula, 'n')
            errorMsg = [errorMsg met ' contains n in formula; '];
        end
        model = addMetabolite(model, [met '[' chosenComp ']'], res.name, res.formula, res.chebi, met, res.pubchem);
        model.metLIPIDMAPSID{end} = res.lipidmaps;
        model.metCharges(end) = 0; % TODO maybe there is a better way
    end

    % put rxn into model
    model=addReaction(model,KEGG_ID,'reactionName',r_name,'reactionFormula',rxnEquation);
    model.rxnEC(end)={strjoin(ecs, ',')};
    addedRxnId = KEGG_ID;
    if ~exist('rxnNote','var')
        rxnNote='';
    end
    model.rxnNotes(end)={rxnNote};
    if ~exist('gpr','var')
        gpr='';
    end
    model.grRules(end)={gpr};
    if ~exist('subSystems','var')
        subSystems=char.empty;
    end
    model.subSystems(end) = {subSystems};
end

% this function always puts a metabolite into [h], but maches metabolites
% to their existing counterpart in [h]
function [rxnEq, newMets, chosenComp, errorMsg] = modifyForModel(rxnEq, model)
    % Find all metabolites
    mets = regexp(rxnEq, '[C][0-9]+', 'match');
    newMets = cell(length(mets), 1);
    chosenComp = 'h';
    errorMsg = '';

    % search in models and find most likely compartment, and possible
    % existing mets
    poss_exist = struct();
    for i = 1:length(mets)
        met = mets{i};
        for j = 1:length(model.mets)
            model_met = model.mets{j};
            model_met_kegg = model.metKEGGID{j};
            if (startsWith(model_met, met) || ...
                    contains(model_met_kegg, met)) && ...
                    endsWith(model_met, '[h]')

                % store the possible corresponding met
                if isfield(poss_exist, met)
                    poss_exist.(met) = [poss_exist.(met); {model_met}];
                else
                    poss_exist.(met) = {model_met};
                end
            end
        end
        if ~isfield(poss_exist, met)
            newMets{i} = met;
            poss_exist.(met) = {[met '[' chosenComp ']']};
        end
    end
    newMets(cellfun(@isempty, newMets)) = [];
    
    % Modify all metabolites
    multipleExistFor = cell(length(mets), 1);
    for i = 1:length(mets)
        currPossExist = poss_exist.(mets{i});

        % check if we have multiple possible matches
        if length(currPossExist) > 1
            disp(currPossExist)
            multipleExistFor{i} = mets{i};
        end

        % replace all metabolites with existing model metabolites
        rxnEq = strrep(rxnEq, mets{i}, currPossExist{1});
    end

    multipleExistFor(cellfun(@isempty, multipleExistFor)) = [];
    if ~isempty(multipleExistFor)
        errorMsg = ['Multiple existing mets for ' ...
            strjoin(multipleExistFor, ', ') '; '];
    end
end

% this function finds the best compartment, but does not match the
% metabolites if the existing metabolite's id is not KEGG_id[comp]
function [rxnEq, newMets, chosenComp] = modifyForModel_old(rxnEq, model, KEGG_ID)
    % Find all metabolites
    mets = regexp(rxnEq, '[C][0-9]+', 'match');
    newMets = cell(length(mets), 1);

    % search in models and find most likely compartment, and possible
    % existing mets
    comp_counts = struct();
    for i = 1:length(mets)
        foundInModel = false;
        met = mets{i};
        for j = 1:length(model.mets)
            model_met = model.mets{j};
            model_met_kegg = model.metKEGGID{j};
            if startsWith(model_met, met) || contains(model_met_kegg, met)
                foundInModel = true;
                % Extract the compartment from model_met
                [match, tokens] = regexp(model_met, '\[(.*?)\]', 'match', 'tokens');
                if ~isempty(match)
                    comp = tokens{1}{1};
                    % Add 1 to the counter for the compartment
                    if isfield(comp_counts, comp)
                        comp_counts.(comp) = comp_counts.(comp) + 1;
                    else
                        comp_counts.(comp) = 1;
                    end
                end
            end
        end
        if ~foundInModel
            newMets{i} = met;
        end
    end
    newMets(cellfun(@isempty, newMets)) = [];

    % Get the compartment with the highest number of occurrences
    allComps = fields(comp_counts);
    if ~isempty(allComps)
        [maxcount, maxcompindex] = max(struct2array(comp_counts));
        if isfield(comp_counts, 'h') && comp_counts.h == maxcount
            chosenComp = 'h';
        else
            chosenComp = allComps{maxcompindex};
        end
    else
        chosenComp = 'h';
    end
    
    % this case should almost never happen
    if ~strcmp(chosenComp, 'h')
        disp(['Encountered non-h best compartment:' chosenComp])
        disp([KEGG_ID ': ' rxnEq])
        disp(comp_counts)
    end

    % Modify all metabolites
    for i = 1:length(mets)
        rxnEq = strrep(rxnEq, mets{i}, [mets{i} '[' chosenComp ']']);
    end
end

function sorted = remGeneralECs(ecs)
    % If a more specific EC number is present, remove a general one

    sorted = sort(ecs);
    remove = false(length(ecs), 1);

    % Iterate over the sorted EC numbers
    for i = length(sorted):-1:1
        currEC = sorted{i};

        % truncate everything after first -
        knownPart = currEC(cummin(currEC ~= '-'));
        
        % sum > 1 because the ec itself will always be true
        remove(i) = contains(currEC, '-') && sum(startsWith(sorted, knownPart)) > 1;
    end

    sorted(remove) = [];
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