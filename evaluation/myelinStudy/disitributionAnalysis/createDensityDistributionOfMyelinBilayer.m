clc; clear; close all;

addpath(genpath(sprintf('..%s..%s..%slibrary',filesep,filesep,filesep)));
addpath(genpath(sprintf('..%s..%s..%stxtFiles' ...
    ,createFilesepStringArray(3))));
constants = readConstantsFile('constants.txt');

resultsDir = sprintf("..%s..%s..%sRESULTS%swholeMyelin" ...
    + "_densityDistributions%s",createFilesepStringArray(5));
if ~exist(resultsDir,'dir')
    warning('The results directory does not exist but will be created.');
    mkdir(resultsDir)
end

directory = "C:\Users\maxoe\Documents\Gromacs\testDatasets\";
folderName = "MYELIN";
layerForm = "\Bilayer\";
fileName = "20221222_MYELIN_TIP4_Bilayer_50water_" ...
    + "shortDurationForDistribution_allAtoms_prd_dt40ps_simTime17_322ns";
groFileName = "step6.5_equilibration";

avogadro = constants.avogadroConstant;
massHydrogenKG = constants.atomicWeightHydrogen/1000/avogadro;
massNitrogenKG = constants.atomicWeightNitrogen/1000/avogadro;
massOxygenKG = constants.atomicWeightOxygen/1000/avogadro;
massCarbonKG = constants.atomicWeightCarbon/1000/avogadro;
massPhosphorusKG = constants.atomicWeightPhosphorus/1000/avogadro;
massSulfurKg = constants.atomicWeightSulfur/1000/avogadro;


atomNamesToIgnore = ["POT" "CLA" "M"];
moleculesToIgnore = ["POT" "CLA"];

numberOfLocations = 500;

%% --------- find out which index belongs to which atom and molecule
fprintf("Determining which index belongs to which kind of atom. \n");
groFileID = fopen(sprintf("%s%s.gro",directory + folderName ...
    + layerForm,groFileName));
groFileLine = fgetl(groFileID);
splittedLineString = strsplit(groFileLine," ");
if strcmp(splittedLineString{1}, "Title")
    groFileLine = fgetl(groFileID);
    atomsCountFromGroFile = str2num(groFileLine); %#ok<ST2NM>
    atomNameArray = strings(1,atomsCountFromGroFile);
    moleculeNameArray = strings(1,atomsCountFromGroFile);
    groFileLine = fgetl(groFileID);
end

atomNameCollection = strings();
moleculeNameCollection = strings();
oldMoleculeIndex = 1;
oldMoleculeName = string.empty();
atomWeightCollection = [];

moleculeAtomCollection = [];
moleculeStartIndex = 1;
moleculeCounter = zeros(1,11);

while groFileLine ~= -1
    splittedLineString = strsplit(groFileLine," ");
    groFileLine = fgetl(groFileID);
    
    splittedLineString = splittedLineString(2:end);
    
    if length(splittedLineString) < 8
        moleculeCounter(moleculeNameCollection == moleculeNameArray( ...
                atomIndex - 1)) = moleculeCounter( ...
                moleculeNameCollection == moleculeNameArray( ...
                atomIndex - 1)) + 1;
        break;
    end
    
    moleculeName = splittedLineString{1};
    
    % check for different line styles
    if length(splittedLineString) == 9
        atomName = splittedLineString{2};
        atomIndex = str2num(splittedLineString{3}); %#ok<ST2NM>
    elseif length(splittedLineString) == 8 && atomIndex >= 9999
        atomAndIndex = splittedLineString{2};
        atomName = atomAndIndex(1:end-5);
        atomIndex = str2double(atomAndIndex(end-4:end));
        if atomIndex > atomsCountFromGroFile
            error(['There is a higher index than there are atoms' ...
                ' in the sample according to gro file.']);
        end
    else
        continue;
    end 
    
    % check for different atom descriptions
    if contains(atomName,'POT') || contains(atomName,'CLA')
        atomNameArray(atomIndex) = atomName(1:3);
    else
        atomNameArray(atomIndex) = atomName(1);
    end
    
    moleculeIndex = getMoleculeIndexFromMoleculeName(moleculeName);
    
    % molecule name collection    
    if contains(moleculeName,"CER24") 
        moleculeName = "GAL";
    elseif contains(moleculeName,"BGAL")
        moleculeName = "BGAL";
    elseif contains(moleculeName,"CHL")
        moleculeName = "CHL";
    elseif contains(moleculeName,"DSPE")
        moleculeName = "DSPE";
    elseif contains(moleculeName,"DOPC")
        moleculeName = "DOPC";
    elseif contains(moleculeName,"SOPS")
        moleculeName = "SOPS";
    elseif contains(moleculeName,"SSM")
        moleculeName = "SSM";
    elseif contains(moleculeName,"PSPI")
        moleculeName = "PSPI";
    elseif contains(moleculeName,"POT")
        moleculeName = "POT";
    elseif contains(moleculeName,"CLA")
        moleculeName = "CLA";
    elseif contains(moleculeName,"SOL")
        moleculeName = "WATER";
    else
        error("Molecule not found");
    end
    
    moleculeNameArray(atomIndex) = moleculeName;
    
    if isempty(oldMoleculeName)
        oldMoleculeName = moleculeName;
    end
    
    
    if oldMoleculeIndex ~= moleculeIndex
        moleculeLength = atomIndex - moleculeStartIndex;
        
        if any([1 3] == moleculeIndex) % Then, the molecule = GAL + BGAL
            if moleculeLength == 151
                moleculeNameArray(moleculeStartIndex:atomIndex-1) ...
                    = "GALS";
            else
                moleculeNameArray(moleculeStartIndex:atomIndex-1) ...
                    = "GALC";
            end
        end
        
        % find GALS and GALC
        if moleculeName ~= "BGAL"
            if ~contains(moleculeNameCollection,moleculeNameArray( ...
                    atomIndex - 1))
                if moleculeNameCollection(1) == ""
                    moleculeNameCollection(1) = moleculeNameArray( ...
                        atomIndex - 1);
                else
                    moleculeNameCollection(end+1) = moleculeNameArray( ...
                        atomIndex - 1); %#ok<SAGROW>
                end
            end
            
            
            moleculeCounter(moleculeNameCollection == moleculeNameArray( ...
                atomIndex - 1)) = moleculeCounter( ...
                moleculeNameCollection == moleculeNameArray( ...
                atomIndex - 1)) + 1;
            
            moleculeStartIndex = atomIndex;
            
        end
        oldMoleculeIndex = moleculeIndex;
    end
    
    
    if ~contains(atomNameCollection,atomNameArray(atomIndex))
        switch atomNameArray(atomIndex)
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
            atomNameCollection(end) = atomNameArray(atomIndex);
        else
            atomNameCollection(end+1) = atomNameArray(atomIndex); %#ok<SAGROW>
        end
        atomWeightCollection(end+1) = weight; %#ok<SAGROW>
    end
    
end

%%
fprintf("Loading data and plotting distributions.\n");
trajData = load(sprintf("%s%s.mat",directory + folderName  ...
    + layerForm,fileName));
densityDistribution = [];
densityDistributionOnlyHdydrogen = [];

[atomsCountFromTraj,dimensions,numberOfTimeSteps] = size( ...
    trajData.trajectories(:,:,1:end));

if atomsCountFromTraj ~= atomsCountFromGroFile
    error(['The number of atoms in the gro file and' ...
        ' the trajectory file mismatch.']);
end

positionsInNM = double(trajData.trajectories(:,:,1:end));
positionsAroundZero = positionsInNM - mean(positionsInNM,1);
positionsInM = positionsAroundZero * constants.nanoMeter;
dXForTimeStep = zeros(1,numberOfTimeSteps);
dYForTimeStep = zeros(1,numberOfTimeSteps);
dZForTimeStep = zeros(1,numberOfTimeSteps);
dVolumeForTimeStep = zeros(1,numberOfTimeSteps);


discreteLocations = zeros(numberOfTimeSteps,numberOfLocations);
frequency_timeLocMolAtom = zeros(numberOfTimeSteps,numberOfLocations ...
    ,length(moleculeNameCollection),length(atomNameCollection));


for timeStepNr = 1:numberOfTimeSteps
    
    if mod(timeStepNr,20) == 0
        fprintf("Time step: %i\n", timeStepNr);
    end
    
    minX = min(positionsInM(:,1,timeStepNr));
    maxX = max(positionsInM(:,1,timeStepNr));
    discreteLocationBorders = linspace(minX,maxX,numberOfLocations+1);
    discreteLocations(timeStepNr,:) = discreteLocationBorders(1:end-1) ...
        + diff(discreteLocationBorders)/2;
    discreteLocationBorders(1) = discreteLocationBorders(1) ...
        - constants.nanoMeter;
    discreteLocationBorders(end) = discreteLocationBorders(end) ...
        + constants.nanoMeter;
    
    
    for atomNr = 1:atomsCountFromTraj
        atomName = atomNameArray(atomNr);
        if any(atomNamesToIgnore == atomName)
            continue;
        end
        moleculeName = moleculeNameArray(atomNr);
        if any(moleculesToIgnore == moleculeName)
            continueM
        end
        
        atomXPos = positionsInM(atomNr,1,timeStepNr);
        locationIndex = find(discreteLocationBorders <= atomXPos,1,'last');
        
        frequency_timeLocMolAtom(timeStepNr,locationIndex ...
            ,moleculeNameCollection == moleculeName ...
            ,atomNameCollection == atomName) = ...
            frequency_timeLocMolAtom(timeStepNr,locationIndex ...
            ,moleculeNameCollection == moleculeName ...
            ,atomNameCollection == atomName) + 1;
        
    end
    
end


%% plotting
avgLoc = mean(discreteLocations,1);
waterIndex = moleculeNameCollection == "WATER";
hydrogenIndex = atomNameCollection == "H";
avgFreqWater = sum(sum(squeeze(mean(frequency_timeLocMolAtom( ...
    :,:,waterIndex,:),1)),3),2);
avgFreqMemb = sum(sum(squeeze(mean(frequency_timeLocMolAtom( ...
    :,:,~waterIndex,:),1)),3),2);
avgFreqWaterH = squeeze(mean(frequency_timeLocMolAtom(:,:,waterIndex,hydrogenIndex),1));
avgFreqMembH = squeeze(sum(mean(frequency_timeLocMolAtom(:,:,~waterIndex,hydrogenIndex),1),3));

fig1 = initializeFigure();
title("Frequency all atoms");
initializeSubplot(fig1,2,1,1);
plot(avgLoc,avgFreqWater);
plot(avgLoc,avgFreqMemb);
legend("Water","Membrane");


initializeSubplot(fig1,2,1,2);
title("Frequency only hydrogen atoms");
plot(avgLoc,avgFreqWaterH);
plot(avgLoc,avgFreqMembH);
legend("Water","Membrane");


saveFigureTo(resultsDir,"wholeMyelin",datestr(now,"yyyymmdd") ...
    ,"atomDistribution",true);

%% saving
save(resultsDir + datestr(now,"yyyymmdd") + "_myelinDistributionData" ...
    ,'dVolumeForTimeStep','dXForTimeStep','dYForTimeStep' ...
    ,'dZForTimeStep','discreteLocations','atomNameArray' ...
    ,'moleculeNameArray','atomNameCollection' ...
    ,'moleculeNameCollection','frequency_timeLocMolAtom','constants');



%% functions
function moleculeIndex = getMoleculeIndexFromMoleculeName(moleculeName)
indexAsCell = textscan(moleculeName,'%7d');
moleculeIndex = indexAsCell{1};

end




