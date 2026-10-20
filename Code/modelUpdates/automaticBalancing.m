function model = automaticBalancing(model)
%AUTOMATICBALANCING automatic reaction balancing with protons and water
    [unbalancedRxns, imbalances] = get_unbalanced(model);
    
    for i = 1:length(unbalancedRxns)
        rid = unbalancedRxns{i};
        imba = imbalances{i};
        
        % If only H, O, charge in imbalance, check if ratios match H20, H+
        if isempty(setdiff(fieldnames(imba), {'H', 'O', 'charge'}))
            if isfield(imba, 'H')
                h = imba.H;
            else
                h = 0;
            end
            if isfield(imba, 'O')
                o = imba.O;
            else
                o = 0;
            end
            if isfield(imba, 'charge')
                c = imba.charge;
            else
                c = 0;
            end
    
            % Adjust reactions based on balance conditions
            if h == 2 * o + c
                % Print if imbalance is suspicously large
                % There is one reaction with 4 protons, where balancing is fine
                if max(abs([o, c])) > 4
                    disp(['Large imbalance in reaction ' rid ': ' ...
                        num2str(o) ', ' num2str(c)]);
                    disp('Not changed')
                else
                    rxnIdx = strcmp(model.rxns, rid);
                    waterIdx = strcmp(model.mets, 'C00001[h]');
                    protonIdx = strcmp(model.mets, 'C00080[h]');
                    model.S(waterIdx, rxnIdx) = model.S(waterIdx, rxnIdx) + o;
                    model.S(protonIdx, rxnIdx) = model.S(protonIdx, rxnIdx) - c;
                end
            end
        end
    end
end

