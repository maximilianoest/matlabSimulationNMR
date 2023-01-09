clc; clear all; close all; fclose("all");

%% add directories
addpath(genpath(sprintf('..%s..%slibrary',filesep,filesep)));
addpath(genpath(sprintf('..%s..%stxtFiles',createFilesepStringArray(2))));
constants = readConstantsFile('constants.txt');
NUMBER_OF_heaviest_MOLECULES_IN_ONE_LEAFLET = 50;

atomNameCollection = ["C" "H" "N" "O" "P" "S"];
atomWeightCollectionInKG = [constants.atomicWeightCarbon ...
    constants.atomicWeightHydrogen ...
    constants.atomicWeightNitrogen ...
    constants.atomicWeightOxygen ...
    constants.atomicWeightPhosphorus ...
    constants.atomicWeightSulfur]/constants.avogadroConstant/1000;

if length(atomNameCollection) ~= length(atomWeightCollectionInKG)
    error("Atom names and weights do not coincide.")
end

%% composition based on Schyboll2019
% number of lipids in one of two leaflets -> overall 100 lipid molecules

lipidClasses = ["CHL" "PC" "PE" "GalC" "GalS" ];
frequencyOfEachLipid = [27 23 23 23 4];
atomCompositionsOfLipids = [ ....
    27 46 0 1 0  0; ...
    40 80 1 8 1  0; ...
    39 76 1 8 1  0; ...
    40 77 1 7 0  0; ...
    40 76 1 11 0 1];

fprintf("According to Schyboll2019 - 'Dipolar induced spin lattice ...' " ...
    + "(10.1038/s41598-019-51003-4) in one leaflet: \n");
schyboll2019 = createGromacSimulationConfiguration( ...
    atomCompositionsOfLipids,lipidClasses,atomNameCollection ...
    ,atomWeightCollectionInKG,frequencyOfEachLipid);
clearvars lipidComposition frequencyOfEachLipid atomCompositionsOfLipids


%% composition based on Schyboll2020
% -> number of lipids in one of two leaflets

lipidClasses = ["CHL" "PC" "PE" "GalC" "PG"];
frequencyOfEachLipid = [13 12 11 12 2];
atomCompositionsOfLipids = [ ...
    27 46 0  1 0 0; ...
    40 80 1  8 1 0; ...
    39 76 1  8 1 0; ...
    40 77 1  7 0 0; ...
    38 75 0 10 1 0 ];

fprintf("According to Schyboll2020 - 'Origin of orientation dep...' " ...
    + "(DOI 10.1002/mrm.28277) in one leaflet: \n");
schyboll2020 = createGromacSimulationConfiguration( ...
    atomCompositionsOfLipids,lipidClasses,atomNameCollection ...
    ,atomWeightCollectionInKG,frequencyOfEachLipid);
clearvars lipidComposition frequencyOfEachLipid atomCompositionsOfLipids

%% composition based on DeVries1981
% percent based on total lipid weight -> summed to 100%
% Ganglioside are neglected because they have an influence < 0.5% where
% even the S.E. is larger.

fprintf("According to DeVries1981 - 'Lipid composition of axolemma...'" ...
    + "(PMID: 7240954):\n") 
lipidClasses = ["CHL"  "PC"  "PE" "GalC" "GalS"  "SM"   "PS"  "PI"];
weightInPercent = [0.271 0.105 0.184 0.214   0.08  0.073  0.063  0.010];
deVries1981 = createConfigurationBasedOnHistology(lipidClasses ...
    ,weightInPercent);

%% composition based on Morell1984
fprintf("According to Rasband2012 based on Morell1984 " ...
    + "- 'Myelin structure ...' (DOI: 10.1016/b978-0-12-374947-5.00010-9): \n");
lipidClasses = ["CHL" "GalC" "GalS" "PE" "PC" "SM" "PS" "PI"];
weightInPercent = [0.277 0.227 0.038 0.156 0.112 0.079 0.048 0.006];
weightNormalizedTo1 = weightInPercent/sum(weightInPercent);
Rasband2012 = createConfigurationBasedOnHistology(lipidClasses ...
    ,weightNormalizedTo1);

%% composition based on Fewster1976
% percent dry weight of tissue (inlcuding proteins and so)
fprintf("According to Fewster1976 - 'Lipid composition ...' " ...
    + "(DOI: 10.1007/bf00313273): \n");
lipidClasses = ["CHL" "PC" "PS" "PI" "SM" "GalC" "GalS" "PE"];
weightInPercent = [0.200 0.069 0.084 0.021 0.048 0.120 0.058 0.137];
weightNormalizedTo1 = weightInPercent/sum(weightInPercent);
Fewster1976 = createConfigurationBasedOnHistology(lipidClasses ...
    ,weightNormalizedTo1);

%% composition based on OBrien1965
% precent dry weight (including lipids and proteins) for 4 different 
% persons with different age -> taking the average over these four people
fprintf("According to O'Brien1965 - ' Lipid composition of the normal...'" ...
    + " (PMID: 5865382) \n");
lipidClasses = ["CHL" "PE" "PS" "PC" "SM" "GalC" "GalS"];
weightInPercent = [ ...
    (18.6+21.5+18.6+19.7) ...
    (14.2+11.3+14.2+11.2) ...
    (5.5+4.2+5.5+5.3) ...
    (12.1+9.1+12.2+8.3) ...
    (4.6+4.4+4.6+4.4) ...
    (13.7+19.2+14.0+16.0) ...
    (5.1+3.9+5.1+3.4)]/4;
weightNormalizedTo1 = weightInPercent/sum(weightInPercent);
oBrien1965 = createConfigurationBasedOnHistology(lipidClasses ...
    ,weightNormalizedTo1);
    

%% lipids and their compositions which are not covered by schyboll

% % WIDDER2018: Normal composition To investigate the interaction of the myelin basic protein with the
% myelin-like lipid monolayer, chloroform solutions of cholesterol, PE, PS, PC, SM, and PI
% were prepared in a mole ratio of 44:27:13:11:3:2 (further referred to as 'normal' composition).
% Thus, the lipid mixture consists of 56% amphiphilic phospholipids and sphingomyelin,
% and of 44% cholesterol. This composition corresponds to that of the cytoplasmic leaet of
% the myelin sheath.

% Some lipids which weren't included in the schbyoll2019 and schyboll2020
% model has to be added. Therefore, their weight need to be determined.
% Based on the lipid studies which were made by Groll/Oesterreicher these
% the lipids DPPC, SPG and DOPS are used. The following compositions are
% found. 

% CHL, GalC GalS Ganglioside, PC, PE, PE PL, PS, SM Plasmalogens, PG

% Based on the word document given by the lipid study team from
% israel:
% - PC 
%   Chosen according to document from israel team:
%   phosphatidylcholine (18:2 – 18 carbon tail length with 2 
%   double bonds(most) + 16:0 – 16 carbon tail length with no double bond )
%   according to deVries1981, this kind of lipid only occurs with 0.3%. For
%   better comparison with Shtangel2020 and previous work 
%   -> comparing with CHARMM => PLPC: 
%   -> composition by counting in strucutural formula in word document from
%   team israel = c42h80n1o8p1 and websites
%   counting in .gro document:
groFilePath = "C:\Users\maxoe\Documents\Gromacs\testDatasets\PLPC\prd.gro";
atomCounterCollection = getNumberOfAtomsInLipidBasedOn(groFilePath ...
    ,atomNameCollection);
weightOfPLPC = atomCounterCollection' * atomWeightCollectionInKG;
% - PS: phosphatidylserine (18:1 – 18:1 carbon tail length with ...
%   1 double bonds)
%   -> comparing with CHARMM => 



%% compare weight distribution of all references

referenceArray = {deVries1981,"deVries1981" ...
    ;Rasband2012,"Rasband2012";Fewster1976,"Fewster1976"; ...
    oBrien1965,"oBrien1965"};

referencesSummary = struct();
for refNr = 1:size(referenceArray,1)
    referenceValues = referenceArray{refNr,1};
    referenceName = referenceArray{refNr,2};
    referencesSummary(refNr).referenceName = referenceName;
    for lipidNameNr = 1:length(referenceValues)
        referencesSummary(refNr).(referenceValues(lipidNameNr).name) ...
            = referenceValues(lipidNameNr).weightInPercent;
    end
end
numberOfReferences = length(referencesSummary);
referencesSummary(end+1).referenceName = "AVERAGE";
averageLipidDistribution = [];
lipidClasses = string(fieldnames(referencesSummary))';
lipidClasses = lipidClasses(2:end);

for lipidName = lipidClasses
    allValues = getAllValuesFromAFieldInStructAsArray(referencesSummary( ...
        1:end-1),lipidName);
    if allValues == "" 
        warning("For the lipid %s there were less than %i values given." ...
            ,lipidName,numberOfReferences);
    end
    referencesSummary(end).(lipidName) = mean(str2double(allValues( ...
        allValues ~= "")));
    averageLipidDistribution(end+1) = mean(str2double(allValues( ...
        allValues ~= ""))); %#ok<SAGROW>
end

%% Determine number of lipid molecules for this distribution

fprintf("The summed up compositions values are %.4f.\nFor each lipid: \n" ...
    ,sum(averageLipidDistribution));
for counter = 1:length(averageLipidDistribution)
    fprintf("%7s : %.4f \n",lipidClasses(counter) ...
        ,averageLipidDistribution(counter));
end

%% lipid compositions of different lipids based on a sample gro file
dataDirectory = "C:\Users\maxoe\Documents\Gromacs\testDatasets\myelinModel\allLipidsComposition\";
fileName = "step5_input.gro";
fileID = fopen(dataDirectory+fileName);
firstLine = fgetl(fileID);
secondLine = fgetl(fileID);

lipidsInGroData = lipid.empty;
line = fgetl(fileID);
lipidIDToSkip = -1;
oldLipidID = -1;
while true
   nextLine = fgetl(fileID);
   if nextLine == -1
       break;
   elseif contains(line,["TIP3" "CLA" "POT"])
       continue;
   end
   
   [lipidID,lipidName,atomName] = getAtomAndLipidNameFromGroFileLine(line);
   
   if lipidID ~= oldLipidID
       
       if lipidID == 1
           oneLipid = lipid(lipidName,length(atomNameCollection));
       elseif lipidName ~= "BGAL"
           for otherLipid = lipidsInGroData
               existingLipidName = otherLipid.lipidName;
               existingAtomComposition = otherLipid.atomComposition;
               if (existingLipidName == oneLipid.lipidName) ...
                       .* (existingAtomComposition == oneLipid.atomComposition) %#ok<*BDSCA,*BDLOG>
                   oneLipid = lipid(lipidName,length(atomNameCollection));
                   break;
               end
           end
           if sum(oneLipid.atomComposition) ~= 0
               lipidsInGroData(end+1) = oneLipid; %#ok<SAGROW>
               oneLipid = lipid(lipidName,length(atomNameCollection));
           end
       end
   end
   
   line = nextLine;
   oldLipidID = lipidID;
   oldLipidName = lipidName;
   
   oneLipid.increaseAtomCompositionCounter( ...
       atomNameCollection == atomName);
   
end

%% collect lipid information
heaviestLipid = lipid.empty();
for lipidNr = 1:length(lipidsInGroData)
    if lipidsInGroData(lipidNr).lipidName == "CER24"
        if lipidsInGroData(lipidNr).atomComposition( ...
                atomNameCollection == "S") > 0
            if lipidsInGroData(lipidNr).atomComposition( ...
                    atomNameCollection == "H") == 92
                newLipidName = "GalS1/0";
            else
                newLipidName = "GalS1/1";
            end
        else
            if lipidsInGroData(lipidNr).atomComposition( ...
                    atomNameCollection == "H") == 93
                newLipidName = "GalC1/0";
            else
                newLipidName = "GalC1/1";
            end
        end
        lipidsInGroData(lipidNr).lipidName = newLipidName;
    end
    
    lipidsInGroData(lipidNr).calculateLipidWeight(atomWeightCollectionInKG);
    lipidsInGroData(lipidNr).determineLipidClass(lipidClasses);
    lipidsInGroData(lipidNr).determineRelativeFrequencyInMyelinModel( ...
        lipidClasses,averageLipidDistribution);
    
    fprintf("Class: %4s Lipid: %4s, weight: %4.4d, frequency: %4.4f \n" ...
        ,lipidsInGroData(lipidNr).lipidClass ...
        ,lipidsInGroData(lipidNr).lipidName ...
        ,lipidsInGroData(lipidNr).lipidWeight ...
        ,lipidsInGroData(lipidNr).relativeFrequencyInMyelinModel);
    
    if isempty(heaviestLipid)
        heaviestLipid = lipidsInGroData(lipidNr);
    elseif lipidsInGroData(lipidNr).lipidWeight > heaviestLipid.lipidWeight
       heaviestLipid = lipidsInGroData(lipidNr); 
    end
end

%% determine composition and write into CSV file

dataDirectory = ...
    "C:\Users\maxoe\Google Drive\Promotion\Simulation\evaluation\myelinComposition\";
fileName = "lipidCompositions_EXPERIMENTAL.csv";
fileID = fopen(dataDirectory + fileName);
header = fgetl(fileID);
splittedHeader = string(strsplit(header,";"));
newLine = string(sprintf("%s;numberOfLipidsInModel;moleculeWeight" ...
    + ";frequencyInModel\n",header));
fclose(fileID);

modelWeightBasedOnHeaviestLipid = NUMBER_OF_heaviest_MOLECULES_IN_ONE_LEAFLET ...
    * heaviestLipid.lipidWeight;
fprintf("The overall weight of the lipid model is %4.4d \n" ...
    ,modelWeightBasedOnHeaviestLipid);

fprintf("Model contains the following numbers of lipid" ...
    + " (one of each class):\n")
lipidCountInModel = [];
for lipidNr = 1:length(lipidsInGroData)
   lipidInCollection = lipidsInGroData(lipidNr);
   lipidCountInModel(end+1) = modelWeightBasedOnHeaviestLipid ...
       * lipidInCollection.relativeFrequencyInMyelinModel ...
       / lipidInCollection.lipidWeight; %#ok<SAGROW>
   fprintf("%5s #: %3i (%4.4f)\n",lipidInCollection.lipidName ...
       ,round(lipidCountInModel(end)),lipidCountInModel(end));
   
   
   fileID = fopen(dataDirectory + fileName);
   firstLine = fgetl(fileID);
   nextLine = fgetl(fileID);
   while nextLine ~= -1
       splittedLine = strsplit(nextLine,";");
       lipidName = splittedLine{splittedHeader == "lipidName"};
       if lipidName == "CHL"
           lipidName = "CHL1";
       end
       lipidClass = splittedLine{splittedHeader == "lipidClass"};
       if lipidInCollection.lipidName == lipidName ...
               && lipidInCollection.lipidClass == lipidClass 
           newLine(end+1) = string(sprintf("%s;%i;%.5d;%.4f\n",nextLine ...
               ,round(lipidCountInModel(end)) ...
               ,lipidInCollection.lipidWeight ...
               ,lipidInCollection.relativeFrequencyInMyelinModel)); %#ok<SAGROW>
           
           lipidsInGroData(lipidNr).numberOfMoleculesInModel ...
               = lipidCountInModel(end);
           lipidsInGroData(lipidNr).frequencyOfLipidWithSameClass ...
               = splittedLine{splittedHeader  ...
               == "relativeFrequencyOfTailLength"};
           break;
       end
       nextLine = fgetl(fileID);
   end
   fclose(fileID);
end

[~,sortingIndices] = sort(newLine(2:end));
sortedNewLines = [newLine(1) newLine(sortingIndices+1)];
fileID = fopen(dataDirectory + "modelCompositionWithLipidCounts.csv",'w');
fprintf(fileID,"%s",sortedNewLines);
fclose(fileID);

%% Determine model composition for simulations
finalCompositionStruct = struct();
for lipidNr = 1:length(lipidsInGroData)
    lipidOfInterest = lipidsInGroData(lipidNr);
    
    for otherLipidNr = 1:length(lipidsInGroData)
        lipidToCompare = lipidsInGroData(otherLipidNr);
        if lipidOfInterest.lipidClass == lipidToCompare.lipidClass
            if lipidOfInterest.frequencyOfLipidWithSameClass ...
                >= lipidToCompare.frequencyOfLipidWithSameClass
                finalCompositionStruct.(lipidOfInterest.lipidClass) ...
                    = lipidOfInterest;
            else
                finalCompositionStruct.(lipidOfInterest.lipidClass) ...
                    = lipidToCompare;
            end
        end
        
    end
end

finalCompositionTable = struct.empty;
moleculeNumberSum = 0;
summedFrequency = 0;
overallWeight = 0;
for fieldName = string(fieldnames(finalCompositionStruct))'
    lipidOfInterest = finalCompositionStruct.(fieldName);
    finalCompositionTable(end+1).lipidClass = lipidOfInterest.lipidClass; %#ok<SAGROW>
    finalCompositionTable(end).lipidName = lipidOfInterest.lipidName;
    finalCompositionTable(end).numberOfMolecules ...
        = round(lipidOfInterest.numberOfMoleculesInModel);
    finalCompositionTable(end).moleculeWeight ...
        = lipidOfInterest.lipidWeight;
    finalCompositionTable(end).relativeFrquency ...
        = lipidOfInterest.relativeFrequencyInMyelinModel;
    moleculeNumberSum = moleculeNumberSum ...
        + round(lipidOfInterest.numberOfMoleculesInModel);
    summedFrequency = summedFrequency ...
        + lipidOfInterest.relativeFrequencyInMyelinModel;
    overallWeight = overallWeight ...
        + lipidOfInterest.numberOfMoleculesInModel ...
        * lipidOfInterest.lipidWeight;
end

finalCompositionTable(end+1).lipidClass = "SUM";
finalCompositionTable(end).relativeFrquency = summedFrequency;
finalCompositionTable(end).numberOfMolecules = moleculeNumberSum;
finalCompositionTable(end).moleculeWeight = overallWeight;


%% check if created files are in compliance with the caluclated ones
% check model in C:\Users\maxoe\Google Drive\Promotion\Simulation\RESULTS 
% \myelinModelComposition\20221222_myelinModelComposition_Bilayer 
% \20221222_myelin_bilayer_50water for completeness and compare to the
% results above.
fileID = fopen("C:\Users\maxoe\Google Drive\Promotion\Simulation\RESULTS\" ...
    + "myelinModelComposition\20221222_myelinModelComposition_Bilayer\" ...
    + "charmm-gui-7171069307_UNTOUCHED\gromacs\step5_input.gro");

lipidsInGroDataCheck = getLipidsInModelFromGroFile(fileID ...
    ,atomNameCollection);

for lipidNr = 1 : length(lipidsInGroDataCheck)
    if lipidsInGroDataCheck(lipidNr).lipidName == "CER24"
        if lipidsInGroDataCheck(lipidNr).atomComposition( ...
                atomNameCollection == "S") > 0
            lipidsInGroDataCheck(lipidNr).lipidName = "GalS1/1";
        else
            lipidsInGroDataCheck(lipidNr).lipidName = "GalC1/1";
        end
    end
end


%% functions
function lipidsInGroData = getLipidsInModelFromGroFile(fileID ...
    ,atomNameCollection) 
lipidsInGroData = lipid.empty;

% first and second line of gro file header
fgetl(fileID); 
fgetl(fileID);

% first usable line in gro file
line = fgetl(fileID);
oldLipidID = -1;
while true
   nextLine = fgetl(fileID);
   if nextLine == -1
       lipidsInGroData(end).numberOfMoleculesInModel ...
           = lipidsInGroData(end).numberOfMoleculesInModel + 1;
       break;
   elseif contains(line,["TIP3" "CLA" "POT"])
       continue;
   end
   
   [lipidID,lipidName,atomName] = getAtomAndLipidNameFromGroFileLine(line);
   
   if lipidID ~= oldLipidID
       
       if lipidID == 1
           oneLipid = lipid(lipidName,length(atomNameCollection));
       elseif lipidName ~= "BGAL"
           for otherLipid = lipidsInGroData
               existingLipidName = otherLipid.lipidName;
               existingAtomComposition = otherLipid.atomComposition;
               if (existingLipidName == oneLipid.lipidName) ...
                       .* (existingAtomComposition == oneLipid.atomComposition) %#ok<*BDSCA,*BDLOG>
                   lipidsInGroData(end).numberOfMoleculesInModel ...
                       = lipidsInGroData(end).numberOfMoleculesInModel + 1;
                   oneLipid = lipid(lipidName,length(atomNameCollection));
                   break;
               end
           end
           if sum(oneLipid.atomComposition) ~= 0
               lipidsInGroData(end+1) = oneLipid; %#ok<AGROW>
               oneLipid = lipid(lipidName,length(atomNameCollection));
           end
       end
   end
   
   line = nextLine;
   oldLipidID = lipidID;
   
   oneLipid.increaseAtomCompositionCounter( ...
       atomNameCollection == atomName);
   
end
end

function atomCounterCollection = getNumberOfAtomsInLipidBasedOn(groFilePath ...
    ,atomNameCollection)
groFileID = fopen(groFilePath);
groFileLine = fgetl(groFileID);

atomCounterCollection = zeros(1,length(atomNameCollection));
while groFileLine ~= -1
    splittedLineString = strsplit(groFileLine," ");
    groFileLine = fgetl(groFileID);
    if length(splittedLineString) < 2
        continue;
    elseif splittedLineString{2} == "1PLPC"
        atomName = splittedLineString{3};
        atomName = atomName(1);
        atomCounterCollection(atomNameCollection == atomName) = ...
            atomCounterCollection(atomNameCollection == atomName) + 1;
    else
        break;
    end
end
fclose(groFileID);
end
function lipidFrequency = getLipidFrequency(structure,lipidName)
structureLipids = getAllValuesFromAFieldInStructAsArray(structure,"name");
index = structureLipids == lipidName;
if sum(index) == 0
    lipidFrequency = [];
    return;
end
lipidFrequency = structure(index).weightInPercent;
end
function lipidWeight = getLipidWeight(structure,lipidName)
structureLipids = getAllValuesFromAFieldInStructAsArray(structure,"name");
index = structureLipids == lipidName;
if sum(index) == 0
    lipidWeight = [];
    return;
end
lipidWeight = structure(index).lipidWeight;
end
function array = getAllValuesFromAFieldInStructAsArray(structure,fieldName)
structureLength = length(structure);
array = string();
for index = 1:structureLength
    
    value = structure(index).(fieldName);
    if isempty(value)
        value = "";
    end
    array(end+1) = value; %#ok<AGROW>
    
end
array = array(2:end);


end
function informationStruct = createConfigurationBasedOnHistology( ...
    lipidComposition,weightInPercent)
if length(lipidComposition) ~= length(weightInPercent)
    error("Information do not fit together.")
end

[lipidComposition, sortingIndex] = sort(lipidComposition);
weightInPercent = weightInPercent(sortingIndex);

informationStruct = struct();
fprintf("Each lipid in weight perecent of total lipid mass: \n");
for lipidNr = 1:length(lipidComposition)
    fprintf("%10s: %.3f \n",lipidComposition(lipidNr) ...
        ,weightInPercent(lipidNr));
    informationStruct(lipidNr).name = lipidComposition(lipidNr);
    informationStruct(lipidNr).weightInPercent = weightInPercent(lipidNr);
end

fprintf(" \n \n");

end
function configurationStruct = createGromacSimulationConfiguration( ...
    atomCompositionsOfLipids,lipidComposition,atomNameCollection ...
    ,atomWeightCollectionInKG,frequencyOfEachLipid)

if diff(size(atomCompositionsOfLipids) ...
        == [length(lipidComposition) length(atomNameCollection)])
    error("The given data on lipids and atoms are not consistent.");
end

[lipidComposition, sortingIndex] = sort(lipidComposition);
atomCompositionsOfLipids = atomCompositionsOfLipids(sortingIndex,:);
frequencyOfEachLipid = frequencyOfEachLipid(sortingIndex);

overallLeafletWeight = 0;
totalNumberOfLipids = sum(frequencyOfEachLipid);
configurationStruct = struct();
for lipidNr = 1:length(lipidComposition)
    configurationStruct(lipidNr).name = lipidComposition(lipidNr);
    configurationStruct(lipidNr).frequency = frequencyOfEachLipid(lipidNr);
    configurationStruct(lipidNr).relFrequency = ...
        frequencyOfEachLipid(lipidNr)/totalNumberOfLipids;
    configurationStruct(lipidNr).lipidComposition = ...
        atomCompositionsOfLipids(lipidNr,:);
    configurationStruct(lipidNr).lipidWeight = ...
        atomCompositionsOfLipids(lipidNr,:) ...
        * atomWeightCollectionInKG';
    configurationStruct(lipidNr).allLipidsOfOneTypeWeight = ...
        configurationStruct(lipidNr).lipidWeight * ...
        frequencyOfEachLipid(lipidNr);
    overallLeafletWeight = overallLeafletWeight ...
        + configurationStruct(lipidNr).allLipidsOfOneTypeWeight;
end

fprintf("  There are %i lipids in the sample. \n" ...
    + "  The overall weight of the leaflet is %.4d kg. \n" ...
    + "  Each lipid in weight percent of total lipid mass " ...
    + "(# lipids in system, rel. frequ. of lipid in system): \n" ...
    ,totalNumberOfLipids,overallLeafletWeight);
for lipidNr = 1:length(lipidComposition)
    configurationStruct(lipidNr).weightInPercent = ...
        configurationStruct(lipidNr).allLipidsOfOneTypeWeight ...
        /overallLeafletWeight;
    fprintf("%10s: %.3f (%2i, %.3f) \n",configurationStruct(lipidNr).name ...
        ,configurationStruct(lipidNr).weightInPercent ...
        ,configurationStruct(lipidNr).frequency ...
        ,configurationStruct(lipidNr).relFrequency );
end

fprintf(" \n \n");

end

