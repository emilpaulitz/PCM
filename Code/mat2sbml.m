clearvars -except gurobiAvailable projDir; clc;
model_version = '1'; %'2';
org = ''; %'Nicotiana_tabacum.'; %'Arabidopsis_thaliana.';

% load model
if isempty(org)  % union model
    load(strcat(projDir, 'Data/pcm/', org, 'pcm.v', ...
            model_version, '.mat'));
else
    load(strcat(projDir, 'Data/pcm/species/', org, ...
        'pcm.v', model_version, '.mat'));
end
pcm = model;

if isfield(pcm, 'A')
    pcm = rmfield(pcm, 'A');
end

if ~isfield(pcm, 'rxnECNumbers') && isfield(pcm, 'rxnEC')
    pcm.rxnECNumbers = strrep(pcm.rxnEC, ',', ';');
end

% Change format of the subsystems array
% This will give warnings with writeCbModel, but will 
% map subsystems correctly in the .xml
for i = 1:length(pcm.rxns)
    s = pcm.subSystemNames(logical(pcm.rxn2subSystem(i, :)));
    pcm.subSystems{i} = s;
end

% write model
if isempty(org)  % union model
    writeCbModel(pcm, 'format','sbml', 'filename', ...
        [projDir, 'Data/pcm/', org, 'pcm.v', model_version]);
else
    writeCbModel(pcm, 'format','sbml', 'filename', ...
        [projDir, 'Data/pcm/species/', org, 'pcm.v', model_version]);
end