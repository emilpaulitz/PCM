clearvars -except wrep_carved carved gurobiAvailable projDir; clc;

importVersion = 1;
modelPath = [projDir 'Data/pcm/species/'];

%% Arabidopsis
load([modelPath, 'Arabidopsis_thaliana.pcm.v' num2str(importVersion) '.mat'], 'model');
pcm = model;

% translate the gene identifiers from NCBI to TAIR IDs, given as locus_tag
% of the gene in the fasta
fastaFile = [projDir ...
    'Data/sequences/Arabidopsis_thaliana.fasta'];
fid = fopen(fastaFile, 'r');

idToLocusMap = containers.Map;

% This replaces the ncbi IDs like lcl|NC_003070.9_prot_NP_001322048.1_6859
% with the TAIR gene names like At5g64300
while ~feof(fid)
    line = fgetl(fid);
    
    if startsWith(line, '>')
        % Extract the gene identifier
        parts = strsplit(line, ' ');
        geneId = parts{1}(2:end); % remove '>' from the start
        % during speciesSpecModels, | was replaced by _
        geneId = strrep(geneId, '|', '_');
        
        % Extract the locus_tag
        locusTagMatch = regexp(line, '\[locus_tag=([^\]]+)\]', 'tokens');
        if ~isempty(locusTagMatch)
            locusTag = locusTagMatch{1}{1};

            % catch old chloroplast names
            switch lower(locusTag)
            case 'arthcp002'
              locusTag='AtCg00020';
            case 'arthcp005'
              locusTag='AtCg00070';
            case 'arthcp006'
              locusTag='AtCg00080';
            case 'arthcp007'
              locusTag='AtCg00120';
            case 'arthcp008'
              locusTag='AtCg00130';
            case 'arthcp009'
              locusTag='AtCg00140';
            case 'arthcp010'
              locusTag='AtCg00150';
            case 'arthcp016'
              locusTag='AtCg00220';
            case 'arthcp017'
              locusTag='AtCg00270';
            case 'arthcp018'
              locusTag='AtCg00280';
            case 'arthcp019'
              locusTag='AtCg00300';
            case 'arthcp021'
              locusTag='AtCg00340';
            case 'arthcp022'
              locusTag='AtCg00350';
            case 'arthcp025'
              locusTag='AtCg00420';
            case 'arthcp026'
              locusTag='AtCg00430';
            case 'arthcp027'
              locusTag='AtCg00440';
            case 'arthcp028'
              locusTag='AtCg00470';
            case 'arthcp029'
              locusTag='AtCg00480';
            case 'arthcp030'
              locusTag='AtCg00490';
            case 'arthcp031'
              locusTag='AtCg00500';
            case 'arthcp032'
              locusTag='AtCg00510';
            case 'arthcp036'
              locusTag='AtCg00550';
            case 'arthcp037'
              locusTag='AtCg00560';
            case 'arthcp038'
              locusTag='AtCg00570';
            case 'arthcp039'
              locusTag='AtCg00580';
            case 'arthcp042'
              locusTag='AtCg00630';
            case 'arthcp049'
              locusTag='AtCg00680';
            case 'arthcp050'
              locusTag='AtCg00690';
            case 'arthcp051'
              locusTag='AtCg00700';
            case 'arthcp052'
              locusTag='AtCg00710';
            case 'arthcp068'
              locusTag='AtCg00890';
            case 'arthcp071'
              locusTag='AtCg01010';
            case 'arthcp074'
              locusTag='AtCg01050';
            case 'arthcp075'
              locusTag='AtCg01060';
            case 'arthcp076'
              locusTag='AtCg01070';
            case 'arthcp077'
              locusTag='AtCg01080';
            case 'arthcp078'
              locusTag='AtCg01090';
            case 'arthcp079'
              locusTag='AtCg01100';
            case 'arthcp080'
              locusTag='AtCg01110';
            case 'arthcp086'
              locusTag='AtCg01250';
            end
            
            % Store the mapping of geneIdentifier to locusTag
            idToLocusMap(geneId) = locusTag;
        end
    end
end
fclose(fid);

% replace gene names in the model and capitalize them
for i = 1:length(pcm.grRules)
    words = strsplit(pcm.grRules{i});
    newWords = {};
    anyChanged = false;
    for j = 1:length(words)
        if ~isempty(words{j})
            currGene = words{j};
            if isempty(words{j})
                continue
            end

            % remove brackets
            numBracketsFront = 0;
            while startsWith(currGene, '(')
                numBracketsFront = numBracketsFront + 1;
                currGene = currGene(2:end);
            end
            numBracketsBack = 0;
            while endsWith(currGene, '(')
                numBracketsBack = numBracketsBack + 1;
                currGene = currGene(1:end-1);
            end

            % make the replacement
            if isKey(idToLocusMap, currGene)
                currGene = idToLocusMap(currGene);
            end
            % necessary for gpr simplification
            currGene = strrep(currGene, '.', '_');


            % re-add brackets
            currGene = [repmat('(', 1, numBracketsFront) currGene];
            currGene = [currGene, repmat(')', 1, numBracketsBack)];

            % catch old chloroplast gene names and problematic characters
            if startsWith(lower(currGene), 'arthcp') || ...
                contains(currGene, '|') || contains(currGene, '&') || ...
                contains(currGene, 'NC_003070.9_prot_NP_001322048.1_6859')
                disp(currGene)
            end

            anyChanged = true;
            newWords{end+1} = capitalizeGene(currGene);
        end
    end
    if anyChanged
        pcm.grRules{i} = makeGprConformGECKO(strjoin(newWords, ' '));
    end
end
pcm.grRules = cellfun(@simplifyGpr, pcm.grRules, UniformOutput=false);

% some manual curation of Arabidopsis genes
% these are all the same, and need to be the same for GECKO to recognize
for i = 1:length(pcm.grRules)
    pcm.grRules{i} = replace(pcm.grRules{i}, 'At1g25083', 'At1g24909');
    pcm.grRules{i} = replace(pcm.grRules{i}, 'At1g25155', 'At1g24909');
    pcm.grRules{i} = replace(pcm.grRules{i}, 'At3g61010', 'At3g61000');
end
pcm = updateFromGrRules(pcm);

% At4g26300: annotated cytosolic
pcm = removeGenesFromModel(pcm, {'At4g26300'});
% the gpr of rubisco got incorrect after that because only the gene was
% removed, but not the whole complex
newRule = ['AtCg00490 and At1g67090 or AtCg00490 and At5g38410 or ' ...
    'AtCg00490 and At5g38420 or AtCg00490 and At5g38430'];
pcm.grRules(strcmp(pcm.rxns, 'RBO_h')) = {newRule};
pcm.grRules(strcmp(pcm.rxns, 'RBC_h')) = {newRule};

% changes for pc generation
pcm.eccodes = pcm.rxnEC;
pcm.rxnECNumbers = pcm.rxnEC;
pcm = removeGenesFromModel(pcm, pcm.genes(~startsWith(pcm.genes, 'At')));

pcm = updateFromGrRules(pcm);

% activate phototrophic conditions
media = readtable([projDir 'Data/analysis/simulationMedia.tsv'], ...
    'FileType', 'text', 'Delimiter', '\t');
media.hetero = media.hetero_str;
media = removevars(media, ["hetero_suc", "hetero_str"]);
for rxnIx = 1:length(media.rxn)
    currRxnIx = strcmp(media.rxn{rxnIx}, pcm.rxns);
    if any(currRxnIx)
        pcm.ub(currRxnIx) = media.photo(rxnIx);
    end
end

% write to disk
if isfield(pcm, 'A')
    pcm = rmfield(pcm, 'A');
end
pcm.description = 'Arabidopsis PCM';

writeCbModel(pcm, 'fileName', ...
    [modelPath, 'Arabidopsis_thaliana.pcm.v' num2str(importVersion+1) '.mat']);
writeCbModel(pcm, 'fileName', ...
    [modelPath, 'Arabidopsis_thaliana.pcm.v' num2str(importVersion+1) '.xml']);

%% Nicotiana tabacum
clearvars -except gurobiAvailable projDir importVersion modelPath; clc;

load([modelPath, 'Nicotiana_tabacum.pcm.v' num2str(importVersion) '.mat'], 'model');
pcm = model;

% remove "lcl_" from the beginning of gene names
for i = 1:length(pcm.grRules)
    currGpr = pcm.grRules{i};
    
    if ~isempty(currGpr)
        % Matches 'lcl|N' followed by anything (non-greedy), until a space or end of string
        pattern = 'lcl_(N[^\s]*($|\s))';
        pcm.grRules{i} = regexprep(currGpr, pattern, '$1');
    end

    words = strsplit(pcm.grRules{i});
    for j = 1:length(words)
        if contains(words{j}, '|')
            names = strsplit(words{j}, '|');
            words{j} = names{1};
        end
    end
    pcm.grRules{i} = strjoin(words, ' ');
end
% after we removed the problematic |, simplify to make it easier for later
f = waitbar(0, 'Formatting GPRs');
for i = 1:length(pcm.grRules)
    waitbar(i/length(pcm.grRules), f, 'Formatting GPRs');
    currGpr = strrep(pcm.grRules{i}, '.', '__DOT__');
    currGpr = simplifyGpr(currGpr);
    pcm.grRules{i} = strrep(currGpr, '__DOT__', '.');
end
close(f);
pcm = updateFromGrRules(pcm);

% Some manual curation of genes
% These were withdrawn by NCBI because their gene model predictions were
% not supported anymore
genesToRemove = {'NW_015926895.1_prot_XP_016474838.1_66823', ...
                 'NW_015945826.1_prot_XP_016488181.1_78840', ...
                 'NW_015922619.1_prot_XP_016471368.1_63718', ...
                 'NW_015813851.1_prot_XP_016497714.1_11814'};
% This one is a proper protein but was falsely annotated to PSII (its only
% reaction), so remove it. 
genesToRemove = [genesToRemove, ...
                 'NW_015929945.1_prot_XP_016477130.1_68894'];
pcm = removeGenesFromModel(pcm, genesToRemove);
% This one (accD) has only 66% pid to Arabidopsis, so replace it manually
pcm.grRules = strrep(pcm.grRules, 'AtCg00500', ...
    '(NC_001879.2_prot_NP_054508.1_84035 or NC_001879.2_prot_NP_054508.1_34)');
pcm = updateFromGrRules(pcm);

% Output to file
pcm.description = 'Tobacco PCM';

writeCbModel(pcm, 'fileName', ...
    [modelPath, 'Nicotiana_tabacum.pcm.v' num2str(importVersion+1) '.mat']);
writeCbModel(pcm, 'fileName', ...
    [modelPath, 'Nicotiana_tabacum.pcm.v' num2str(importVersion+1) '.xml']);