function model = updateFromGrRules(model)
    %UPDATEFROMGRRULES update genes, rxnGeneMat, and rules fields
    %   Uses information from grRules field to fill the other gene-related
    %   fields
    regenerateGeneNames = false;
    if isfield(model, 'geneNames')
        prevGenes = model.genes;
        prevGeneNames = model.geneNames;
        model = rmfield(model, 'geneNames');
        regenerateGeneNames = true;
    end
    regenerateUniprot = false;
    if isfield(model, 'geneUniprotID')
        prevGenes = model.genes;
        prevGeneUniprot = model.geneUniprotID;
        model = rmfield(model, 'geneUniprotID');
        regenerateUniprot = true;
    end

    if isfield(model, 'rxnGeneMat')
        model = rmfield(model, 'rxnGeneMat');
    end
    if isfield(model, 'genes')
        model = rmfield(model, 'genes');
    end
    if isfield(model, 'rules')
        model = rmfield(model, 'rules');
    end

    res = warning(); % silence the expected warning
    warning('off','all')
    model = updateGenes(model);
    warning(res)

    model = generateRules(model);
    model = buildRxnGeneMat(model);
    model = updateGenes(model);

    if regenerateGeneNames
        model.geneNames = cell(length(model.genes), 1);
        for i = 1:length(prevGeneNames)
            newIx = strcmp(model.genes, prevGenes{i});
            if any(newIx)
                model.geneNames(newIx) = ...
                    prevGeneNames(i);
            end
        end
        for i = 1:length(model.geneNames)
            if isempty(model.geneNames{i})
                model.geneNames{i} = '';
            end
        end
    end
    if regenerateUniprot
        model.geneUniprotID = cell(length(model.genes), 1);
        for i = 1:length(prevGeneUniprot)
            newIx = strcmp(model.genes, prevGenes{i});
            if any(newIx)
                model.geneUniprotID(newIx) = prevGeneUniprot(i);
            end
        end
        for i = 1:length(model.geneUniprotID)
            if isempty(model.geneUniprotID{i})
                model.geneUniprotID{i} = '';
            end
        end
    end
end