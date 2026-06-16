function s = capitalizeGene(x)
%CAPITALIZEGENE Capitalize Arabidopsis genes to conform to Uniprot
% convention: capital A, capital C/M for the chromosome, everything else
% lowercase

    % remove potential brackets in front of the gene
    numBrackets = 0;
    while startsWith(x, '(')
        numBrackets = numBrackets + 1;
        x = x(2:end);
    end

    % bring into UniProt format (capitalized except for C-chromosome)
    if isempty(x)
        s = x;
    else
        if startsWith(lower(x), 'atcg')
            s = ['AtCg' lower(x(5:end))];
        elseif startsWith(lower(x), 'atmg')
            s = ['AtMg' lower(x(5:end))];
        else
            s = [upper(x(1)) lower(x(2:end))];
        end
    end

    % add back brackets
    for i = 1:numBrackets
        s = ['(' s];
    end
end

