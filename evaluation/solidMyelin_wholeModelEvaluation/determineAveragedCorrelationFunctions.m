clear all; close all; clc; fclose('all');

addpath(genpath(sprintf("..%s..%slibrary",filesep,filesep)));
addpath(genpath(sprintf("..%s..%stxtFiles",createFilesepStringArray(2))));

resultsDir = sprintf( ...
    "..%s..%sRESULTS%ssolidMyelin_20230109_DetermineCorrelationFunctions" ...
    ,createFilesepStringArray(3));
directoryInfo = getOnlyFoldersAndDataIn(resultsDir);
corrFuncDirs = directoryInfo(getValuesFromStructAsArray(directoryInfo ...
    ,"isdir"));
loadCorrFuncZerothOrder = 1;
loadCorrFuncFirstOrder = 1;
loadCorrFuncSecondOrder = 1;
specificDate = "20220401";
atomIndices = 'all';
nearestNeighbours = 5500;
simDurationInNs = 50;
dTInPs = 20;

for directoryNr = 1:length(corrFuncDirs)
    directoryName = string(corrFuncDirs(directoryNr).name);
    folderName = string(corrFuncDirs(directoryNr).folder);
    corrFuncDirectory = folderName + filesep + directoryName + filesep;
    switch directoryName
        case specificDate + "_corrFuncZerothOrder"
            if loadCorrFuncZerothOrder
                fileNamesToLoad = getFileNamesOfInterest( ...
                    corrFuncDirectory,atomIndices,nearestNeighbours ...
                    ,simDurationInNs,dTInPs);
                avgCorrFuncZerothOrder = loadAverageCorrFunc( ...
                    folderName+filesep+directoryName,fileNamesToLoad);
            end
        case specificDate + "_corrFuncFirstOrder"
            if loadCorrFuncFirstOrder
                fileNamesToLoad = getFileNamesOfInterest( ...
                    corrFuncDirectory,atomIndices,nearestNeighbours ...
                    ,simDurationInNs,dTInPs);
                avgCorrFuncFirstOrder = loadAverageCorrFunc( ...
                    folderName+filesep+directoryName,fileNamesToLoad);
            end
        case specificDate +"_corrFuncSecondOrder"
            if loadCorrFuncSecondOrder
                fileNamesToLoad = getFileNamesOfInterest( ...
                    corrFuncDirectory,atomIndices,nearestNeighbours ...
                    ,simDurationInNs,dTInPs);
                avgCorrFuncSecondOrder = loadAverageCorrFunc( ...
                    folderName+filesep+directoryName,fileNamesToLoad);
            end
        otherwise
           error("Not detectable directory.");
   end
    
end

save(sprintf("%s%s%s_%s_%s_%s.mat",resultsDir,filesep,specificDate ...
    ,datestr(now,'yyyymmdd_HHMM'),"averageCorrelationFunctions") ...
    ,"avgCorrFuncZerothOrder", "avgCorrFuncFirstOrder" ...
    ,"avgCorrFuncSecondOrder", "nearestNeighbours", "atomIndices" ...
    ,"simDurationInNs", "dTInPs");



%% functions
% 1
function onlyFoldersAndDataInDir = getOnlyFoldersAndDataIn( ...
    absoluteDirectory)

allInDir = dir(absoluteDirectory);
onlyFoldersAndDataInDir = struct.empty();
for element = allInDir' 
    if string(element.name) == "." || string(element.name) == ".."
       continue; 
    end
    if isempty(onlyFoldersAndDataInDir)
        onlyFoldersAndDataInDir = element;
        continue;
    end
    onlyFoldersAndDataInDir(end+1) = element; %#ok<AGROW>
    
end
end

% 2
function valuesFromStructAsArray = getValuesFromStructAsArray( ...
    struct,fieldName)

if isempty(struct)
    error("Input struct is empty. Returned empty array");
    valuesFromStructAsArray = string.empty();
end 

dataType = class(struct(1).(fieldName));
switch dataType
    case 'string'
        valuesFromStructAsArray = string.empty();
    case 'double'
        valuesFromStructAsArray = double.empty();
    case 'single'
        valuesFromStructAsArray = single.empty();
    case 'logical'
        valuesFromStructAsArray = logical.empty();
    otherwise 
        error("Not supported data type. Fix function.");
end

for elementNr = 1:length(struct)
   valuesFromStructAsArray(elementNr) = struct(elementNr).(fieldName);
end


end

%3
function fileNamesToLoad = getFileNamesOfInterest(directory ...
    ,atomIndices,NN,simDurInNs,deltaTInPs)
allFilesInDirectory = what(directory);
matFileNamesInDir = string(allFilesInDirectory.mat);
if length(NN) > 1 || isempty(NN)
    error("The nearest neighbour parameter habe to be a single integer.");
end

fileNamesToLoad = string.empty();
for fileName = matFileNamesInDir'
    splittedName = strsplit(strrep(fileName,".mat",""),"_");
    
    atomIndexFromFileName = str2double(strrep(splittedName( ...
        contains(splittedName,"atomIndex")),"atomIndex",""));
    throwErrorIfVariabelWasNotDeterminableFromFileName( ...
        atomIndexFromFileName);
    if sum(atomIndices ~= 'all') && ...
            sum(atomIndices == atomIndexFromFileName) == 0
        continue;
    end
    
    NNFromFileName = str2double(strrep(splittedName( ...
        contains(splittedName,"NN")),"NN",""));
    throwErrorIfVariabelWasNotDeterminableFromFileName(NNFromFileName);
    if NN ~= NNFromFileName
        continue;
    end
    
    simDurInNsFromFileName = getTimeValueFromFileName(splittedName ...
        ,"simDur","ns");
    throwErrorIfVariabelWasNotDeterminableFromFileName(simDurInNs);
    if simDurInNs ~= simDurInNsFromFileName
        continue;
    end
    
    dTInPsFromFileName = getTimeValueFromFileName(splittedName,"dT","ps");
    throwErrorIfVariabelWasNotDeterminableFromFileName(dTInPsFromFileName);
    if dTInPsFromFileName ~= deltaTInPs
        continue;
    end
    
    fileNamesToLoad(end+1) = fileName; %#ok<AGROW>
end

if isempty(fileNamesToLoad)
    error("No fitting file names were found.");
end

end

% 4
function throwErrorIfVariabelWasNotDeterminableFromFileName(variable)

if isnan(variable)
    error("%s was not determinable from file name.", inputname(1));
end

end

% 5
function numberFromFileName = getTimeValueFromFileName( ... 
    splittedName,stringBeforeNumber,stringAfterNumber)
stringBeforeNumber = char(stringBeforeNumber);
stringAfterNumber = char(stringAfterNumber);
commaInName = sum(abs(contains(splittedName,stringBeforeNumber) ...
    - contains(splittedName,stringAfterNumber)));
if commaInName
    numberBeforeComma = str2double(replace(splittedName(contains( ...
        splittedName,stringBeforeNumber)),stringBeforeNumber,""));
    numberAfterComma = replace(splittedName(contains(splittedName ...
        ,stringAfterNumber)),stringAfterNumber,"");
    numberFromFileName = numberBeforeComma + str2double( ...
        numberAfterComma)/(10^numel(numberAfterComma));
    
else
    numberFromFileName = str2double(replace(splittedName( ...
        contains(splittedName,stringBeforeNumber)) ...
        ,{stringBeforeNumber stringAfterNumber},""));
end
end

% 6
function averagedCorrFunc = loadAverageCorrFunc(directoryPath ...
    ,fileNamesToLoad)

sumCorrFunc = double.empty;
for fileName = fileNamesToLoad
    filePath = directoryPath + filesep + fileName;
    if ~exist(filePath,'file')
        % maybe as array function
        error("File %s does not exist in directory %s." ...
            ,fileName,directoryPath)
    end
    load(filePath,'corrFunc');
    if isempty(sumCorrFunc)
        sumCorrFunc = corrFunc;
        continue;
    end
    sumCorrFunc = sumCorrFunc + corrFunc;
end
fprintf("Averaged %i data curves. \n",length(fileNamesToLoad));
averagedCorrFunc = corrFunc / length(fileNamesToLoad);
end




