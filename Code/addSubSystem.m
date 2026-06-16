function model = addSubSystem(model, rxnList, subsName)
%ADDSUBSYSTEM This function adds a new subsystem to the model, and adds
%this subsystem to the reactions provided in rxnList.
% rxnList has to be an object of the form {'rxn1', 'rxn2'}.
% subsName is the subsystems' name
nSubs = length(model.subSystemNames);
subsIx = nSubs + 1;
model.rxn2subSystem(:, subsIx) = zeros(1, length(model.rxns));
model.subSystemNames{subsIx} = subsName;
for i=1:length(rxnList)
    currRxn = find(strcmp(model.rxns, rxnList{i}));
    if isempty(currRxn)
        disp(['Warning: reaction ', rxnList{i}, ' not found']);
    end
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
    model.rxn2subSystem(strcmp(model.rxns, rxnList{i}), subsIx) = 1;
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