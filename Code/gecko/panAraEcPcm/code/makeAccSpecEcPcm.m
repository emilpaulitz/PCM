% Generate multiple accession-specific ecGems simultaneously
% This script requests input of the root password at some point to be able
% to move the DLKcat.tsv file
clearvars -except gurobiAvailable projDir; clc;

generateNew = false;
runDLKcatNew = true;
checkInstalls = false;
defaultKCat = 26.5;

%% check installations
if checkInstalls
    checkInstallation;
    GECKOInstaller.install;
end

%% Start
dataPath = [projDir 'Code/gecko/panAraEcPcm/data/'];
inModelPath = [projDir 'Data/analysis/accSpecPCM/04_models/'];
outModelPath = [projDir 'Code/gecko/panAraEcPcm/models/'];

fileList = dir(fullfile(inModelPath, '*.mat'));
media = readtable([projDir 'Data/analysis/simulationMedia.tsv'], ...
    'FileType', 'text', 'Delimiter', '\t');

problems = '';
for accIx = 1:length(fileList)
    fileName = fileList(accIx).name;
    acc = erase(fileName, '.mat');
    accUnderscore = strrep(acc, '-', '_');

    if ~generateNew && exist([outModelPath 'Ath.' accUnderscore '.ecModel.tuned.xml'], 'file')
        disp(['Skipping ' acc ' because its model exists'])
        continue
    end

    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    disp(['Starting accession ' acc '!'])
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    
    %% load .mat, curate, export to .xml
    % load acc-spec .mat model
    load([inModelPath, acc, '.mat'], 'model');
    
    model.metNames(strcmp(model.mets, 'C00111[h]')) = {'glycerone phosphate'};

    % curations following what GECKO found for the Ath model from KEGG
    model.eccodes(strcmp(model.rxns, 'R13563')) = {'5.2.1.14'};
    model.eccodes(contains(model.rxns, 'ABCAT_h')) = {'7.-.-.-'};

    % remove all genes not from arabidopsis or this accession
    toRemove = model.genes(~startsWith(lower(model.genes), 'at'));
    model = removeGenesFromModel(model, toRemove);

    % capitalize arabidopsis genes and change logical structure so that
    % GECKO can understand it
    model.grRules = cellfun(@(s) strrep(s, '.', '_'), model.grRules, 'UniformOutput', false);
    for i = 1:length(model.grRules)

        % capitalize arabidopsis genes
        allWords = strsplit(model.grRules{i});
        for wordIx = 1:length(allWords)
            word = allWords{wordIx};
            if (length(word) == 9) && startsWith(lower(word), 'at') && ...
                (lower(word(4)) == 'g') && (lower(word(3)) ~= 'c')
                capitalized = lower(word);
                capitalized = ['A', capitalized(2:end)];
                model.grRules{i} = strrep(model.grRules{i}, word, capitalized);
            end
        end

        % long GPRs cannot be handled by matlab. Leave them for now
        if sum(model.rxnGeneMat(i, :)) > 50
            disp(['Skip formatting of rxn ' model.rxns{i} ' (Idx: ' ...
                num2str(i) ') with ' num2str(sum(model.rxnGeneMat(i, :))) ...
                ' genes'])
            continue
        end

        model.grRules{i} = makeGprConformGECKO(model.grRules{i});
    end
    model = updateFromGrRules(model);

    % make GECKO find the EC numbers
    model.rxnECNumbers = model.eccodes; 

    % activate phototrophic conditions
    for rxnIx = 1:length(media.rxn)
        currRxnIx = strcmp(media.rxn{rxnIx}, model.rxns);
        if any(currRxnIx)
            model.ub(currRxnIx) = media.photo(rxnIx);
        end
    end
    
    % write pcm to disk as sbml, and perform necessary format changes
    for i = 1:length(model.rxns)
        s = model.subSystemNames(logical(model.rxn2subSystem(i, :)));
        model.subSystems{i} = s;
    end
    if isfield(model, 'A')
        model = rmfield(model, 'A');
    end
    writeCbModel(model, 'fileName', ...
        [projDir 'Code/gecko/panAraEcPcm/models/' accUnderscore '.pcm.v2.xml']);
    
    %% set the appropriate adapter file:
    adapterLocation = fullfile(projDir, ['Code/gecko/panAraEcPcm/adapters/', ...
        accUnderscore '/panAraEcPcmAdapter.m']);
    ModelAdapterManager.setDefault(adapterLocation);
    
    %% load the conventional pcm 
    model = loadConventionalGEM();
    
    %% create draft full ecGEM
    fname = [dataPath 'uniprot_' acc '.tsv'];
    % put the file into uniprot.tsv
    copyfile(fname, [dataPath 'uniprot.tsv'])

    % create draft ecModel
    [ecModel, noUniprot] = makeEcModel(model);
    if ~isempty(noUniprot)
        problem = [acc ': noUniprot has ' num2str(length(noUniprot)) ...
            ' entries (first five): ' ...
            strjoin(noUniprot(1:min(5, length(noUniprot))), ', ')];
        disp(problem)
        problems = [problems problem '\n'];
    end

    %% gather eccodes
    ecModel.eccodes = cellfun(@(s) regexprep(strrep(strrep(s, ' ', ''), ...
        ',', ';'), '^;', ''), ecModel.eccodes, 'UniformOutput', false);
    ecModel = getECfromGEM(ecModel);
    noEC = cellfun(@isempty, ecModel.ec.eccodes);
    % this might not be possible; check what was done for Col-0 and repeat
    % manually/ integrate into the GEM
    ecModel = getECfromDatabase(ecModel, noEC);
    
    %% query BRENDA
    kcatList_fuzzy = fuzzyKcatMatching(ecModel);
    
    %% DLKcat
    [ecModel, noSMILES] = findMetSmiles(ecModel);
    
    % writeDLKcatInput
    % the input file belongs to root for some reason
    if exist([dataPath 'DLKcat.tsv'], 'file')
        system(sprintf('sudo chmod %s %s', '666', [dataPath 'DLKcat.tsv']));
    end
    if exist([dataPath acc '_DLKcat.tsv'], 'file') && ~runDLKcatNew
        copyfile([dataPath acc '_DLKcat.tsv'], [dataPath 'DLKcat.tsv']);
    else
        % these are all the default arguments, except for overwrite = true
        writeDLKcatInput(ecModel, true(numel(ecModel.ec.rxns),1), ...
            ModelAdapterManager.getDefault(), ...
            true, ...
            fullfile(dataPath,'DLKcat.tsv'), ...
            true);
        % if you get a Cannot write to destination error here, change the
        % line with movefile in runDLKcat.m to 
        % movefile(fullfile(params.path,'/data/tempDLKcatOutput.tsv'), filePath, 'f');
        runDLKcat();
        system(sprintf('sudo chmod %s %s', '666', [dataPath 'DLKcat.tsv']));
        copyfile([dataPath 'DLKcat.tsv'], [dataPath acc '_DLKcat.tsv']);
    end
    kcatList_DLKcat = readDLKcatOutput(ecModel);
    
    %% Merge BRENDA and DLKcat outputs
    kcatList_merged = mergeDLKcatAndFuzzyKcats(kcatList_DLKcat, ...
        kcatList_fuzzy);
    
    % and assign to ecModel
    ecModel = selectKcatValue(ecModel, kcatList_merged);
    
    % use average kcats for isozymes (otherwise, missing kcats will allow very
    % high rates)
    ecModel = getKcatAcrossIsozymes(ecModel);

    saveEcModel(ecModel, ['ath.' accUnderscore '.tmpVer.ecModel.yml']);
    ecModel = loadEcModel(['ath.' accUnderscore '.tmpVer.ecModel.yml']);

    % use average kcat and MW for missing enzymes
    [ecModel, rxnsMissingGPR, standardMW, standardKcat] = ...
        getStandardKcat(ecModel);
    % actually, use the same default kcat for all accessions
    ecModel.ec.kcat(ismember(ecModel.ec.rxns, rxnsMissingGPR)) = defaultKCat;
    
    % apply to ecModel
    ecModel = applyKcatConstraints(ecModel);
    
    % apply custom kcats as data/customKcats.tsv; check my curation in AraPcm
    % and add any predictions obtained from DLKcat and other tools TODO
    sqd1RxnIx = strcmp(ecModel.ec.rxns, 'SQD1');
    sqd1Enzs = ecModel.ec.enzymes(logical(ecModel.ec.rxnEnzMat(sqd1RxnIx, :)));
    if length(sqd1Enzs) ~= 1
        problem = [acc ': sqd1 ' num2str(length(sqd1Enzs)) ' Enzymes'];
        disp(problem)
        problems = [problems problem '\n'];
    else
        generateCustomKcatFile(acc, sqd1Enzs{1}, dataPath);
        [ecModel, rxnUpdated, notMatch] = applyCustomKcats(ecModel);
    end
    
    %% set ptot
    ecModel = setProtPoolSize(ecModel);
    
    %% Tune too low kcats until Rubisco is the most used protein
    ecModel = setProtPoolSize(ecModel);
    [ecModel, tunedKcats] = sensitivityTuning(ecModel);
    %struct2table(tunedKcats)
    
    % Revert all kcats that came after rubisco, including rubisco
    rbcIx = find(startsWith(tunedKcats.rxns, 'RBC_h') | ...
        startsWith(tunedKcats.rxns, 'RBO_h'), 1 );
    for i = rbcIx:length(tunedKcats.rxns)
        currRxnIx = strcmp(ecModel.ec.rxns, tunedKcats.rxns{i});
        ecModel.ec.kcat(currRxnIx) = tunedKcats.oldKcat(i);
    end
    
    % apply new Kcats
    ecModel = applyKcatConstraints(ecModel);
    ecModel = setProtPoolSize(ecModel);
    
    %% Write to disk
    saveEcModel(ecModel, ['ath.' accUnderscore '.ecModel.yml']);
    ecModel = loadEcModel(['ath.' accUnderscore '.ecModel.yml']);
    
    ecModel.id = ['Ath ' acc 'ecPCM'];
    exportModel(ecModel, [outModelPath 'Ath.' accUnderscore '.ecModel.tuned.xml']);

    %% clean up
    % prevent problem when reading new adapters
    rmpath([projDir, 'Code/gecko/panAraEcPcm/adapters/', accUnderscore '/']);
end

function generateCustomKcatFile(acc, newEnz, dataPath)
    % copy the template to an accession-specific version
    copyfile([dataPath 'customKcats.template.tsv'],...
        [dataPath acc '_customKcats.tsv'])

    % read the freshly copied file
    opts = delimitedTextImportOptions('Delimiter', '\t',...
        'VariableTypes', 'string', 'VariableNamesLine', 1);
    kcatTable = readtable([dataPath acc '_customKcats.tsv'], opts);

    % mofidy the protein
    rowIdx = find(strcmp(kcatTable.rxns, 'SQD1'));
    kcatTable.proteins(rowIdx) = newEnz;

    % write new table where GECKO will find it
    writetable(kcatTable, [dataPath 'customKcats.tsv'], ...
        'Delimiter', '\t', 'FileType', 'text', 'WriteVariableNames', false);
end