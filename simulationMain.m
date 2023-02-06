%% Initialize system
clear all; close all; clc;

addpath(genpath('library'));
addpath(genpath('txtFiles'));
addpath(genpath('scripts'));
addpath(genpath('matFiles'));
constantsFileName = 'constants.txt';
configuration = readConfigurationFile('configMain.txt');
matlabSimulationDate = datestr(date,'yyyymmdd');
[dataDirectory,resultsFileSavingPath,logFilePath] = ...
    setUpDependenciesBasedOn(configuration);
resultsDirectory = getDirectoryOfAbsoluteFilePath(resultsFileSavingPath);

logMessage('System initialized. Dependendencies were added.',logFilePath)

%% Analyse data set based on file name

analyseFileNameAndCreateVariablesInBaseWorkspace(configuration ...
    ,logFilePath);
if whichLipid == "MYELIN"
   dataDirectory = [dataDirectory formOfLayer filesep];
   if ~isfolder(dataDirectory)
       error("Data directory not found.");
   end
end

%% Define constants
checkIfFileExists(constantsFileName,logFilePath);
logMessage('Defining constants.',logFilePath,true);
constants = readConstantsFile('constants.txt');

dipolDipolConstant = 3/4*(constants.vaccumPermeability/(4*pi) ...
    *constants.hbar*constants.gyromagneticRatioOfHydrogenAtom^2)^2 ...
    /(constants.nanoMeter^6);
omega0 = constants.gyromagneticRatioOfHydrogenAtom ...
    *configuration.mainMagneticField;

%% Define simulation parameters
dataFilePath = sprintf('%s%s.mat',dataDirectory,fileName);
checkIfFileExists(dataFilePath,logFilePath);

logMessage('Defining simulation parameters.',logFilePath,true);
[numberOfHs,~,timeSteps] = size(matfile(dataFilePath),'trajectories');
randomSequenceOfAtoms = randperm(numberOfHs);
logMessage(sprintf(['    Found %d hydrogen atoms at %d time ' ...
    'steps of %.3d s'],numberOfHs,timeSteps,deltaTInS),logFilePath,false);

atomsToCalculate = configuration.atomsToCalculate;
logMessage(sprintf('    Will calculate %i atoms',atomsToCalculate) ...
    ,logFilePath,false);

averagingRegionForSpectralDensity = ...
    getValuesFromStringEnumeration( ...
    configuration.averagingRegionForSpectralDensity,';','numeric');

%% Try prallocation first and then load data
% If implemented, the first part of the scriptFileToRun is executed. In the
% try block the scriptFileToRun is executed. If the arrays are preallocated
% the called scriptFileToRun should quit after preallocaiton and then the
% data will be loaded. If this functionality is not implemented, data
% is load first and then the script is executed.

try
    run(configuration.scriptFileToRun);
    logMessage(sprintf(['The script %s has an implementation for' ...
        ' preallocation of arrays first and then load dataset.'] ...
        ,configuration.scriptFileToRun),logFilePath);
catch ME
    if strcmp(ME.identifier,'MATLAB:UndefinedFunction') ...
            && contains(ME.message,'trajectoryX')
        logMessage(sprintf(['The script %s has no implementation for' ...
            ' preallocate arrays first and then load dataset.'] ...
            ,configuration.scriptFileToRun),logFilePath);
    else
        rethrow(ME)
    end
end

if ~configuration.dataLoaded
    logMessage('Start loading data.',logFilePath);
    [trajectoryX,trajectoryY,trajectoryZ,gromacsSimulationConfiguration] ...
        = loadTrajectoriesAndSimConfig(dataFilePath);
    logMessage('Loading data finished.',logFilePath);
else
    logMessage('Data was already loaded in the run before.',dataFilePath);
end
run(configuration.scriptFileToRun);






