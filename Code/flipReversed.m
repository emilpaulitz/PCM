function model = flipReversed(model, rid)
    rxnIx = strcmp(model.rxns, rid);
    prev_lb = model.lb(rxnIx);
    prev_ub = model.ub(rxnIx);
    model.lb(rxnIx) = -prev_ub;
    model.ub(rxnIx) = -prev_lb;
    model.S(:, rxnIx) = - model.S(:, rxnIx);
end

