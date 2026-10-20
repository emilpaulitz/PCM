%% Convert all .mat models to sbml for MEMOTE analysis
clearvars -except gurobiAvailable projDir; clc;

folderPath = ['/home/emil/Desktop/PhD-Synch/' ...
    'pan_chloroplast/pan_chl_model/PCM/Data/pcm/species'];
filePattern = fullfile(folderPath, '*.mat');
files = dir(filePattern);

for k = 1:length(files)
    baseFileName = files(k).name;
    fullFileName = fullfile(folderPath, baseFileName);

    load(fullFileName);
    
    % ATTENTION DO WE GET model OR pcm?
    pcm = model;
    
    if isfield(pcm, 'A')
        pcm = rmfield(pcm, 'A');
    end
    
    if ~isfield(pcm, 'rxnECNumbers') && isfield(pcm, 'rxnEC')
        pcm.rxnECNumbers = strrep(pcm.rxnEC, ',', ';');
    end
    
    
    % Change format of the subsystems array
    % This will give warnings with writeCbModel, but will 
    % map subsystems correctly in the .xml and for plotting purposes
    for i = 1:length(pcm.rxns)
        s = pcm.subSystemNames(logical(pcm.rxn2subSystem(i, :)));
        pcm.subSystems{i} = s;
    end
    
    % write model
    writeCbModel(pcm, 'format','sbml', 'filename', ...
        strrep(fullFileName, '.mat', '.xml'));
end
