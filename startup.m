%% User-specific paths
% If you have access to gurobi, add the correct path to gurobi, otherwise
% leave it empty. For example:
% gurobiPath = '/opt/gurobi1103/linux64/matlab';
gurobiPath = '';
% Change the path to cobratoolbox installation
pathCobraToolbox = '';

%% global variables
% projDir
projDir = pwd;
projDir = [projDir '/'];
gurobiAvailable = isempty(gurobiPath);

%% setup
% cobratoolbox
if isempty(pathCobraToolbox)
    error(['cobratoolbox is required for this project. Please set the ' ...
          'appropriate path in startup.m!']);
end
cd(pathCobraToolbox);
initCobraToolbox
cd(projDir)

% gurobi
if gurobiAvailable
    addpath gurobiPath
    gurobi_setup
end

% paths for this project
addpath("Code/")
addpath("Code/panGenomeAnalysis/")
addpath("Code/modelComparison/")
addpath("Code/validation/")
addpath("Code/gecko/ecAraPcm/")
addpath("Code/gecko/panAraEcPcm/")
addpath("Code/gecko/panAraEcPcm/code/")
addpath("Code/accSpecPCM/")
addpath("Code/cyanoCompare/")