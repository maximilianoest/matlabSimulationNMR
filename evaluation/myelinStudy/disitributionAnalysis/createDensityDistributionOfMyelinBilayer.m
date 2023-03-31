clc; clear all; close all;

addpath(genpath(sprintf('..%s..%s..%slibrary',filesep,filesep,filesep)));
addpath(genpath(sprintf('..%s..%s..%stxtFiles',createFilesepStringArray(3))));
constants = readConstantsFile('constants.txt');

resultsDir = sprintf("..%s..%s..%sRESULTS%swholeMyelin" ...
    + "_densityDistributions",createFilesepStringArray(4));
if ~exist(resultsDir,'dir')
    warning('The results directory does not exist but will be created.');
    mkdir(resultsDir)
end

directory = "C:\Users\maxoe\Documents\Gromacs\testDatasets\";
folderName = "MYELIN";
layerForm = "\Bilayer\";
fileName = "20221222_MYELIN_TIP4_Bilayer_50water" ...
    + "_shortDurationForDistribution_prd_dt2ps_simDur17_322ns";

groFileName = "step5_input";

avogadro = constants.avogadroConstant;
massHydrogenKG = constants.atomicWeightHydrogen/1000/avogadro;
massNitrogenKG = constants.atomicWeightNitrogen/1000/avogadro;
massOxygenKG = constants.atomicWeightOxygen/1000/avogadro;
massCarbonKG = constants.atomicWeightCarbon/1000/avogadro;
massPhosphorusKG = constants.atomicWeightPhosphorus/1000/avogadro;
massSulfurKg = constants.atomicWeightSulfur/1000/avogadro;


atomNamesToIgnore = ["POT" "CLA" "M"];

numberOfLocations = 500;

%% --------- find out which index belongs to which atom and molecule
fprintf("Determining which index belongs to which kind of atom. \n");
groFileID = fopen(sprintf("%s%s.gro",directory + folderName ...
    + layerForm,groFileName));
groFileLine = fgetl(groFileID);
atomNameCollection = strings();
atomWeightCollection = [];
while groFileLine ~= -1
    splittedLineString = strsplit(groFileLine," ");
    groFileLine = fgetl(groFileID);
    
    % check for title
    if strcmp(splittedLineString{1}, "Generated")
        atomsCountFromGroFile = str2num(groFileLine); %#ok<ST2NM>
        atomNameArray = strings(1,atomsCountFromGroFile);
        moleculeIndices = strings(1,atomsCountFromGroFile);
        groFileLine = fgetl(groFileID);
        continue;
    end
    splittedLineString = splittedLineString(2:end);
    
    molecule = splittedLineString{1};
    % check for different line styles
    if length(splittedLineString) == 6
        atomName = splittedLineString{2};
        index = str2num(splittedLineString{3}); %#ok<ST2NM>
    elseif length(splittedLineString) == 5 && index >= 9999
        atomAndIndex = splittedLineString{2};
        atomName = atomAndIndex(1:end-5);
        index = str2double(atomAndIndex(end-4:end));
        if index > atomsCountFromGroFile
            error(['There is a higher index than there are atoms' ...
                ' in the sample according to gro file.']);
        end
    else
        continue;
    end 
    
    moleculeIndices(index) = molecule;
    % check for different atom descriptions
    if contains(atomName,'POT') || contains(atomName,'CLA')
        atomNameArray(index) = atomName(1:3);
    else
        atomNameArray(index) = atomName(1);
    end
    
    % check for collection of atoms and their weight
    if ~contains(atomNameCollection,atomNameArray(index))
        switch atomNameArray(index)
            case "H"
                weight = massHydrogenKG;
            case "O"
                weight = massOxygenKG;
            case "P"
                weight = massPhosphorusKG;
            case "C"
                weight = massCarbonKG;
            case "N"
                weight = massNitrogenKG;
            case "S"
                weight = massSulfurKg;
            case cellstr(atomNamesToIgnore)
                weight = 1000;
            otherwise
                error("Unknown atom name");
        end
        if atomNameCollection(1) == ""
            atomNameCollection(end) = atomNameArray(index);
        else
            atomNameCollection(end+1) = atomNameArray(index); %#ok<SAGROW>
        end
        atomWeightCollection(end+1) = weight; %#ok<SAGROW>
    end
    
end


fprintf("Loading data and plotting distributions.\n");
trajData = load(sprintf("%s%s.mat",directory + folderName  ...
    + layerForm,fileName));
densityDistribution = [];
densityDistributionOnlyHdydrogen = [];

[atomsCountFromTraj,dimensions,numberOfTimeSteps] = size(trajData.trajectories);

if atomsCountFromTraj ~= atomsCountFromGroFile
    error(['The number of atoms in the gro file and' ...
        ' the trajectory file mismatch.']);
end

positionsInNM = double(trajData.trajectories);
positionsAroundZero = positionsInNM - mean(positionsInNM,1);
positionsInM = positionsAroundZero * constants.nanoMeter;
dXForTimeStep = zeros(1,numberOfTimeSteps);
dYForTimeStep = zeros(1,numberOfTimeSteps);
dZForTimeStep = zeros(1,numberOfTimeSteps);
dVolumeForTimeStep = zeros(1,numberOfTimeSteps);
pos_atomCounter = zeros(size(dXForTimeStep,2) ...
    ,numberOfLocations,length(moleculesToCompare) ...
    ,length(atomNameCollection));
locations = zeros(dimensions,numberOfLocations);

for timeStepNr = 1:size(positionsInNM,3)
    if mod(timeStepNr,200) == 0
        fprintf("Time step: %i\n", timeStepNr);
    end
end




