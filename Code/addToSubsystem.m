function model = addToSubsystem(model, rxnList, subsName)
    if ~any(strcmp(model.subSystemNames, subsName))
        disp(['Did not find the specified subsystem' subsName])
    end
    for i = 1:length(rxnList)
        currRxn = find(strcmp(model.rxns, rxnList{i}));

        if cellIsEmpty(model.subSystems{currRxn})
            model.subSystems{currRxn} = subsName;
        else
            if iscell(model.subSystems{currRxn})
                prev = model.subSystems{currRxn};
            else
                prev = model.subSystems(currRxn);
            end
            model.subSystems(currRxn) = strcat(prev, ',', subsName);
        end
    
        model.rxn2subSystem(currRxn, ...
            strcmp(model.subSystemNames, subsName)) = 1;
    end
end

function isEmpty = cellIsEmpty(cellEle)
    if isempty(cellEle)
        isEmpty = true;
    elseif iscell(cellEle)
        % If it's a cell, check its contents recursively
        isEmpty = all(cellfun(@cellIsEmpty, cellEle));
    else
        isEmpty = false;
    end
end