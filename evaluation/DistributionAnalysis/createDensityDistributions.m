clc; clear all; close all; fclose('all');

% ATTENTION: with mean positions, it is possible that more than one atom is
% at one position. Thus the distribution and density is not right.
% Therefore, a snapshot of the production run is used.

% Beacuse the positions of the lipid molecules are above the boundaries of
% the system when "whole" is chosen, noJump is the better choice. Not sure
% if prd.trr data are joJump but it seems like that.

addpath(genpath(sprintf('..%s..%slibrary',filesep,filesep)));
addpath(genpath(sprintf('..%s..%stxtFiles',createFilesepStringArray(2))));
constants = readConstantsFile('constants.txt');

resultsDir = sprintf( ...
    '..%s..%sRESULTS%sdensityDistributions' ...
    ,createFilesepStringArray(3));
if ~exist(resultsDir,'dir')
    error('The results directory does not exist.');
end
addpath(genpath(resultsDir));

savingDirectory = sprintf('%s%sPlots%s',resultsDir ...
    ,createFilesepStringArray(2));
saving = 1;

%% get relevant results from files
directory = "C:\Users\maxoe\Documents\Gromacs\testDatasets\";
lipidNames = ["DOPS" "PLPC" "PSM"];
filesToLoad = [ ...
    "20220110_DOPS_TIP4_Monolayer_50water_prd_dt10ps_simTime10ns" ...
    "20220804_PLPC_TIP4_Monolayer_50water_prd_dt10ps_simTime10ns" ...
    "20220804_PSM_TIP4_Monolayer_50water_prd_dt10ps_simTime10ns"];

groFileName = "prd";

% fig = initializeFigure();
legendEntries = {};
avogadro = constants.avogadroConstant;
massHydrogenKG = constants.atomicWeightHydrogen/1000/avogadro;
massNitrogenKG = constants.atomicWeightNitrogen/1000/avogadro;
massOxygenKG = constants.atomicWeightOxygen/1000/avogadro;
massCarbonKG = constants.atomicWeightCarbon/1000/avogadro;
massPhosphorusKG = constants.atomicWeightPhosphorus/1000/avogadro;

atomNamesToIgnore = ["POT" "CLA" "M"];

numberOfLocations = 100;

for datasetNr = 1:length(lipidNames)
    lipidName = lipidNames(datasetNr);
    moleculesToCompare = [lipidName, "SOL"];
%% --------- find out which index belongs to which atom and molecule
    disp("Determining which index belongs to which kind of atom");
    groFileID = fopen(sprintf("%s%s%s%s.gro",directory,lipidName...
        ,filesep,groFileName));
    groFileLine = fgetl(groFileID);
    atomNameCollection = strings();
    atomWeightCollection = [];
    while groFileLine ~= -1
       splittedLineString = strsplit(groFileLine," "); 
       groFileLine = fgetl(groFileID);
       % check for title
       if strcmp(splittedLineString,"Title")
           atomsCountFromGroFile = str2num(groFileLine); %#ok<ST2NM>
           atomIndices = strings(1,atomsCountFromGroFile);
           moleculeIndices = strings(1,atomsCountFromGroFile);
           groFileLine = fgetl(groFileID);
           continue;
       end
       splittedLineString = splittedLineString(2:end);
       
       molecule = splittedLineString{1};
       % check for different line styles
       if length(splittedLineString) == 9
           atomName = splittedLineString{2};
           index = str2num(splittedLineString{3}); %#ok<ST2NM>
       elseif length(splittedLineString) == 8 && index >= 9999
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
           atomIndices(index) = atomName(1:3);
       else
           atomIndices(index) = atomName(1);
       end
       
       % check for collection of atoms and their weight
       if ~contains(atomNameCollection,atomIndices(index))
           switch atomIndices(index)
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
               case cellstr(atomNamesToIgnore)
                   weight = 1000;
               otherwise 
                   error("Unknown atom name");
           end
           if atomNameCollection(1) == ""
                atomNameCollection(end) = atomIndices(index);
           else
               atomNameCollection(end+1) = atomIndices(index); %#ok<SAGROW>
           end
           atomWeightCollection(end+1) = weight; %#ok<SAGROW>
       end
       
    end

    
%% --------- load trajectory data and plotting distributions
    disp("Loading data and plotting distributions");
    trajData = load(sprintf("%s%s%s%s.mat",directory,lipidName ...
        ,filesep,filesToLoad(datasetNr))); 
    densityDistribution = [];
    densityDistributionOnlyHdydrogen = [];

    [atomsCountFromTraj,~,timeSteps] = size(trajData.trajectories);
    
    if atomsCountFromTraj ~= atomsCountFromGroFile
        error(['The number of atoms in the gro file and' ...
            ' the trajectory file mismatch.']);
    end
    
    positionsInNM = double(trajData.trajectories);
    positionsAroundZero = positionsInNM - mean(positionsInNM,1);
    positionsInM = positionsAroundZero * constants.nanoMeter;
    dXForTimeStep = zeros(1,size(positionsInNM,3));
    dYForTimeStep = zeros(size(dXForTimeStep));
    dZForTimeStep = zeros(size(dXForTimeStep));
    dVolumeForTimeStep = zeros(size(dXForTimeStep));
    pos_atomCounter = zeros(size(dXForTimeStep,2) ...
        ,numberOfLocations,length(moleculesToCompare) ...
        ,length(atomNameCollection));
    locations = zeros(size(dXForTimeStep,2),numberOfLocations);
    
    for timeStepNr = 1:size(positionsInNM,3)
        if mod(timeStepNr,100) == 0
            fprintf("Time step: %i\n", timeStepNr);
        end
        locations(timeStepNr,:) = linspace(min( ...
            positionsInM(:,1,timeStepNr)),max(positionsInM( ...
            :,1,timeStepNr)),numberOfLocations);
        pos_molName_atomName_Counter = zeros(numberOfLocations ...
            ,length(moleculesToCompare),length(atomNameCollection));
        for atomNr = 1:atomsCountFromTraj
            atomName = atomIndices(atomNr);
            moleculeName = moleculeIndices(atomNr);

            if contains(atomName,atomNamesToIgnore)
               continue; 
            end

            atomXPos = positionsInM(atomNr,1,timeStepNr);
            atomYPos = positionsInM(atomNr,2,timeStepNr);
            atomZPos = positionsInM(atomNr,3,timeStepNr);

            for locationNr = 1:size(locations,2)-1
                if locationNr == 1
                    lowerBoundDecrease = - constants.nanoMeter;
                    upperBoundIncrease = 0;
                elseif locationNr == size(locations,2)-1
                    lowerBoundDecrease = 0;
                    upperBoundIncrease = constants.nanoMeter;
                else
                    lowerBoundDecrease = 0;
                    upperBoundIncrease = 0;
                end
                lowerBound = locations(timeStepNr,locationNr)  ...
                    + lowerBoundDecrease;
                upperBound = locations(timeStepNr,locationNr ...
                    + 1) + upperBoundIncrease;
                if atomXPos >= lowerBound && atomXPos < upperBound
                    % find the right index for the molecule type
                    for moleculeNameNr = 1:length(moleculesToCompare)
                       if contains(moleculeName ...
                               ,moleculesToCompare(moleculeNameNr))
                           moleculeNameCollectionIndex = moleculeNameNr;
                           break;
                       end
                    end
                    % find the right atom name
                    atomNameCollectionIndex = find( ...
                        atomNameCollection == atomName);
                    if length(atomNameCollectionIndex) > 1
                        error('too many');
                    end

                    % sum up the number of atoms within a specific location
                    pos_molName_atomName_Counter(locationNr ...
                        ,moleculeNameCollectionIndex ...
                        ,atomNameCollectionIndex) = ...
                        pos_molName_atomName_Counter(locationNr ...
                        ,moleculeNameCollectionIndex ...
                        ,atomNameCollectionIndex) + 1;
                    continue;
                end
            end  
        end

        if atomsCountFromTraj < sum(sum(sum(pos_molName_atomName_Counter)))
            error("Too many!");
        end

        dX = mean(diff(locations(timeStepNr,:)));
        dY = max(positionsInM(:,2,timeStepNr)) ...
            - min(positionsInM(:,2,timeStepNr));
        dZ = max(positionsInM(:,3,timeStepNr)) ...
            - min(positionsInM(:,3,timeStepNr));
        dVolume = dX*dY*dZ;
        
        dXForTimeStep(timeStepNr) = dX; 
        dYForTimeStep(timeStepNr) = dY; 
        dZForTimeStep(timeStepNr) = dZ; 
        dVolumeForTimeStep(timeStepNr) = dVolume;
        
        pos_atomCounter(timeStepNr,:,:,:) = pos_molName_atomName_Counter;
        
        for moleculeOfInterestNr = 1:length(moleculesToCompare)
            atomDependentLocationsForEachMolecule = ...
                squeeze(pos_molName_atomName_Counter( ...
                :,moleculeOfInterestNr,:));
            massDistribution = atomDependentLocationsForEachMolecule ...
                * atomWeightCollection';
            massDistributionOnlyHydrogen = ...
                atomDependentLocationsForEachMolecule ...
                * (atomWeightCollection.*(atomNameCollection == "H"))';
            densityDistribution(timeStepNr,moleculeOfInterestNr,:) ...
                = massDistribution' / dVolume; %#ok<SAGROW>
            densityDistributionOnlyHdydrogen(timeStepNr ...
                ,moleculeOfInterestNr,:) = ...
                massDistributionOnlyHydrogen'/dVolume; %#ok<SAGROW>
        end
    end
    
    initializeFigure();
    avgLocations = mean(locations,1);
    avgDensityDistribution = squeeze(mean(densityDistribution,1));
    legendEntries = {};
    for moleculeOfInterestNr = 1:length(moleculesToCompare)
        densityDistributionForMolecule = avgDensityDistribution( ...
            moleculeOfInterestNr,:);        
        plot(avgLocations,densityDistributionForMolecule);
        if moleculesToCompare(moleculeOfInterestNr) == "SOL"
            legendEntries{end+1} = "Water"; %#ok<SAGROW>
        else
            legendEntries{end+1} = moleculesToCompare( ...
                moleculeOfInterestNr); %#ok<SAGROW>
        end
    end
    
    title("Density distribution for all atoms in system");
    xlabel('Position $[nm]$');
    ylabel('Density $[kg/m^3]$');
    lgd = legend(legendEntries);
    lgd.Location = "east";
    
    %%
    if saving
        saveFigureTo(savingDirectory,"LocationDistributions",lipidName...
            ,"wholeSystem");
        save(sprintf("%s%s%sdensityDistribution_%s" ...
            ,savingDirectory,filesep,datestr(now,"yyyymmdd"),lipidName) ...
            ,"avgDensityDistribution","avgLocations" ...
            ,"densityDistribution","densityDistributionOnlyHdydrogen" ...
            ,"atomWeightCollection","atomNameCollection" ...
            ,"atomsCountFromGroFile","moleculesToCompare" ...
            ,"dXForTimeStep","dYForTimeStep","dZForTimeStep" ...
            ,"dVolumeForTimeStep","locations","lipidName" ...
            ,"pos_atomCounter");
    end
    
    initializeFigure();
    avgLocations = mean(locations,1);
    avgDensityDistributionOnlyForHydrogen = squeeze(mean(...
        densityDistributionOnlyHdydrogen,1));
    legendEntries = {};
    for moleculeOfInterestNr = 1:length(moleculesToCompare)
        densityDistributionForMolecule = ...
            avgDensityDistributionOnlyForHydrogen(moleculeOfInterestNr,:);        
        plot(avgLocations,densityDistributionForMolecule);
        if moleculesToCompare(moleculeOfInterestNr) == "SOL"
            legendEntries{end+1} = "Water"; %#ok<SAGROW>
        else
            legendEntries{end+1} = moleculesToCompare( ...
                moleculeOfInterestNr); %#ok<SAGROW>
        end
    end
    
    title("Density distribution only for hydrogen atoms");
    xlabel('Position $[nm]$');
    ylabel('Density $[kg/m^3]$');
    lgd = legend(legendEntries);
    lgd.Location = "east";
    
    if saving
        saveFigureTo(savingDirectory,"LocationDistributions",lipidName...
            ,"onlyHydrogen");
    end
    
end
     



