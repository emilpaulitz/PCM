function model = mergeMets(model, obsoleteMetId, keptMetId, annoFieldsList)
    obsIx = strcmp(model.mets, obsoleteMetId);
    keptIx = strcmp(model.mets, keptMetId);

    % merge metabolites in all rxns
    rxnsWithBoth = (model.S(keptIx, :) ~= 0) & (model.S(obsIx, :));
    if any(rxnsWithBoth)
        rxnIx = find(rxnsWithBoth);
        disp(['Potential problem merging ' obsoleteMetId ' and ' ...
            keptMetId ': They occur in the same rxns, this is the first Ix:' ...
            num2str(rxnIx(1))]);
    end
    model.S(keptIx, :) = model.S(keptIx, :) + model.S(obsIx, :);

    % merge annotations in case they exist and are not equal
    for i = 1:length(annoFieldsList)
        annoField = annoFieldsList{i};
        obsAnno = model.(annoField)(obsIx);
        keptAnno = model.(annoField)(keptIx);
        % numeric elements like charges cannot be concatenated, so just
        % keep the kept entry
        if ~(isnumeric(obsAnno) || isnumeric(keptAnno)) && ...
            ~isempty(obsAnno{1}) && ~strcmp(obsAnno{1}, keptAnno{1})
                model.(annoField)(keptIx) = {strcat(keptAnno{1}, ',', ...
                    obsAnno{1})};
        end
    end

    % delete obsolete met
    [model, ~] = removeMetabolites(model, obsoleteMetId, false);
end