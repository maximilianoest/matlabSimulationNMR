logMessage(sprintf('Starting the external script: %s.m \n\n', mfilename) ...
    ,logFilePath,true);


%% Preallocation
if ~exist('nearestNeighbourCases','var')
    logMessage('Preallocation of arrays as a first part of the script' ...
        ,logFilePath,true);
    nearestNeighbourCases = getValuesFromStringEnumeration( ...
        configuration.nearestNeighbourCases,";","numeric");
    nearestNeighbourCases = sort(nearestNeighbourCases,'descend');
    maxNearestNeighbours = nearestNeighbourCases(1);
    
    averagingRegionForSpecDens = getAveragingRegionForSpecDens( ...
        configuration.averagingRegionForSpectralDensity,logFilePath);
    
    atomCounter = 1;
    calculatedAtomIndices = zeros(1,atomsToCalculate);
    
    randomSequenceOfAtoms = randperm(numberOfHs);
    
    nearestNeighboursIDs = zeros(1,maxNearestNeighbours);
    nearestNeighbourDistancesPow3 = zeros(maxNearestNeighbours,timeSteps ...
        ,'single');
    
    theta = configuration.fibreAnglesTheta;
    yAxis = [0 1 0];
    
    phi = configuration.fibreAnglesPhi;
    zAxis = [0 0 1];
    
    polarAngle = zeros(maxNearestNeighbours,timeSteps,'single');
    azimuthAngle = zeros(maxNearestNeighbours,timeSteps,'single');
    
    sphericalHarmonicZerothOrder = zeros(maxNearestNeighbours,timeSteps ...
        ,'like',double(1j));
    sphericalHarmonicFirstOrder = zeros(maxNearestNeighbours,timeSteps ...
        ,'like',double(1j));
    sphericalHarmonicSecondOrder = zeros(maxNearestNeighbours,timeSteps ...
        ,'like',double(1j));
    
    r1Estimation = zeros(1,numberOfHs);
    
    corrFuncZerothOrder = zeros(1,timeSteps,'like',double(1j));
    corrFuncFirstOrder = zeros(1,timeSteps,'like',double(1j));
    corrFuncSecondOrder = zeros(1,timeSteps,'like',double(1j));
    
    sumCorrFuncZerothOrder = zeros(1,timeSteps,'like',double(1j));
    sumCorrFuncFirstOrder = zeros(1,timeSteps,'like',double(1j));
    sumCorrFuncSecondOrder = zeros(1,timeSteps,'like',double(1j));
    
    corrFuncDirectories = createCorrFuncDirectoriesIfNotExist( ...
        resultsDirectory,gromacsSimulationDate);
    logMessage("Created directories for correlation functions" ...
        + " if necessary.",logFilePath,false);
    
    logMessage('Arrays successfully preallocated',logFilePath,false);
    return;
end
meanPositions = single([mean(trajectoryX,2) mean(trajectoryY,2) ...
        mean(trajectoryZ,2)]);
atomTimer = zeros(1,atomsToCalculate);

%% Start simulation
logMessage(sprintf("You have chosen theta = %.4f and phi = %.4f",theta,phi) ...
    ,logFilePath,false);
logMessage("Rotating dataset.",logFilePath);
rotationMatrixPhi = get3DRotationMatrix(deg2rad(phi),zAxis);
rotationMatrixTheta = get3DRotationMatrix(deg2rad(theta),yAxis);
totalRotationMatrix = rotationMatrixTheta*rotationMatrixPhi;
[trajectoryX,trajectoryY,trajectoryZ]  ...
    = rotateTrajectoriesWithRotationMatrix(totalRotationMatrix ...
    ,trajectoryX,trajectoryY,trajectoryZ);

logMessage('Start simulation.',logFilePath,true);


for atomNumber = randomSequenceOfAtoms(atomCounter:atomsToCalculate)
    atomTimerStart = tic;
    calculatedAtomIndices(atomCounter) = atomNumber; %#ok<SAGROW>
    
    [trajectoryX,trajectoryY,trajectoryZ] = calculateRelativePositions( ...
        trajectoryX,trajectoryY,trajectoryZ,atomNumber);
    [nearestNeighboursIDs,nearestNeighbourDistancesPow3] =  ...
        findNearestNeighboursIDs(maxNearestNeighbours,trajectoryX ...
        ,trajectoryY,trajectoryZ);
    
    [polarAngle,azimuthAngle] = ...
        transformToSphericalCoordinates( ...
        trajectoryX(nearestNeighboursIDs,:) ...
        ,trajectoryY(nearestNeighboursIDs,:) ...
        ,trajectoryZ(nearestNeighboursIDs,:));            
    [sphericalHarmonicZerothOrder,sphericalHarmonicFirstOrder ...
        ,sphericalHarmonicSecondOrder] = calculateSphericalHarmonics( ...
        polarAngle,azimuthAngle,nearestNeighbourDistancesPow3);
    
    % calculate correlation functions
    corrFuncZerothOrder = ...
        calculateCorrFuncForMultipleNNCasesAsMatrix( ...
        sphericalHarmonicZerothOrder,nearestNeighbourCases);
    corrFuncFirstOrder = ...
        calculateCorrFuncForMultipleNNCasesAsMatrix( ...
        sphericalHarmonicFirstOrder,nearestNeighbourCases);
    corrFuncSecondOrder = ...
        calculateCorrFuncForMultipleNNCasesAsMatrix( ...
        sphericalHarmonicSecondOrder,nearestNeighbourCases);
    
    % saving
    whichDirectory = resultsDirectory + corrFuncDirectories( ...
        corrFuncDirectories == gromacsSimulationDate + "_" ...
        + getNamesOfVariblesAsArray(corrFuncZerothOrder) + filesep);
    if ~isfolder(whichDirectory)
        error("Chosen wrong directory: %s",whichDirectory);
    end
    saveCorrFuncTo(whichDirectory,corrFuncZerothOrder,atomNumber ...
        ,nearestNeighbourCases,simulationDurationInNs,deltaTInPs ...
        ,theta,phi);
    
    
    whichDirectory = resultsDirectory + corrFuncDirectories( ...
        corrFuncDirectories == gromacsSimulationDate + "_" ...
        + getNamesOfVariblesAsArray(corrFuncFirstOrder) + filesep);
    if ~isfolder(whichDirectory)
        error("Chosen wrong directory: %s",whichDirectory);
    end
    saveCorrFuncTo(whichDirectory,corrFuncFirstOrder,atomNumber ...
        ,nearestNeighbourCases,simulationDurationInNs,deltaTInPs...
        ,theta,phi);
    
    
    whichDirectory = resultsDirectory + corrFuncDirectories( ...
        corrFuncDirectories == gromacsSimulationDate + "_" ...
        + getNamesOfVariblesAsArray(corrFuncSecondOrder) + filesep);
    if ~isfolder(whichDirectory)
        error("Chosen wrong directory: %s",whichDirectory);
    end
    saveCorrFuncTo(whichDirectory,corrFuncSecondOrder,atomNumber ...
        ,nearestNeighbourCases,simulationDurationInNs,deltaTInPs...
        ,theta,phi);
    
    
    sumCorrFuncZerothOrder = sumCorrFuncZerothOrder + corrFuncZerothOrder;
    sumCorrFuncFirstOrder = sumCorrFuncFirstOrder + corrFuncFirstOrder;
    sumCorrFuncSecondOrder = sumCorrFuncSecondOrder + corrFuncSecondOrder;
    
    [averageSpectralDensityFirstOrder,averageSpectralDensitySecondOrder] ...
        = calculateSpectralDensities(sumCorrFuncFirstOrder/atomCounter ...
        ,sumCorrFuncSecondOrder/atomCounter,omega0,deltaTInS ...
        ,averagingRegionForSpecDens);
    r1Estimation(atomCounter) = calculateR1WithSpectralDensity( ...
        averageSpectralDensityFirstOrder ...
        ,averageSpectralDensitySecondOrder,dipolDipolConstant); %#ok<SAGROW>
    
    logMessage(sprintf('Calculated %i atom(s)',atomCounter),logFilePath);
    atomTimer(atomCounter) = toc(atomTimerStart);
    averageTimeForOneAtom = seconds(mean(atomTimer(1:atomCounter)));
    logMessage(sprintf(['  --> Average time for one atom: %s \n' ...
       '        Approximately ready on: %s.'] ...
        ,datestr(averageTimeForOneAtom,'HH:MM:SS') ...
        ,datetime('now')+averageTimeForOneAtom ...
        *(atomsToCalculate-atomCounter)),logFilePath,false);
    
    if mod(atomCounter,configuration.savingIntervall) == 0
        lastSavingDate = datestr(now,'yyyymmdd_HHMM');
        createDataSavingObject();
        save(resultsFileSavingPath,'-struct','dataSavingObject','-v7.3');
        logMessage('Saved data',logFilePath,false);
        
        logMessage(sprintf( ...
            '   => Estimated relaxation rate for calculated atoms: %.4f' ...
            ,r1Estimation(atomCounter)),logFilePath,false);
    end
    
    printEqualSignBreakLineToLogFile(logFilePath);
    atomCounter = atomCounter + 1;
end

