function [unbalancedRxns, imbalances] = get_unbalanced(model)
% A function to get the imbalances for all non-pseudo reactions

    unbalancedRxns = {};
    imbalances = {};
    last_idx = 1;
    
    % Use cobratoolbox checkMassChargeBalance function
    [massImbalance, imBalancedMass, imBalancedCharge, ...
        imBalancedRxnBool, elements, missingFormulaeBool, ...
        balancedMetBool] = checkMassChargeBalance(model, 0);

    % Loop over reactions in the  model
    for i = 1:length(model.rxns)
        reactionID = model.rxns{i};
        
        % Skip pseudo reactions
        if contains(lower(reactionID), 'biomass') || ...
           contains(lower(reactionID), '_im') || ...
           contains(lower(reactionID), 'im_') || ...
           contains(lower(reactionID), '_ex') || ...
           contains(lower(reactionID), 'ex_') || ...
           contains(lower(reactionID), 'exch_') || ...
           contains(lower(reactionID), 'sk_') || ...
           contains(lower(reactionID), 'dm_') || ...
           contains(reactionID, 'Pool') || ...
           contains(lower(reactionID), 'slime')
            continue;
        end
        
        % If reaction was found to be imbalanced
        if imBalancedRxnBool(i)
            % collect reaction IDs of imbalanced reactions
            unbalancedRxns{last_idx} = reactionID;

            % skip reactions flagged because of missing formulas
            currMassImba = imBalancedMass(i);
            if strcmp(currMassImba, 'NaN')
                imbalances{last_idx} = struct();
                last_idx = last_idx + 1;
                continue
            end

            % parse imbalances
            tokens = regexp(currMassImba{1}, ...
                '(-?\d+)\s([A-Z]+)', 'tokens');
            numbers = cellfun(@(x) str2double(x{1}), tokens);
            elements = cellfun(@(x) x{2}, tokens, ...
                'UniformOutput', false);
            curr_imbalances = cell2struct(num2cell(numbers), elements, 2);
            curr_imbalances(1).charge = imBalancedCharge(i);
            imbalances{last_idx} = curr_imbalances;
            last_idx = last_idx + 1;
        end
    end
end