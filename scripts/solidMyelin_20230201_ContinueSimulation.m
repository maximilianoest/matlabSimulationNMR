logMessage(sprintf('Starting the external script: %s.m \n\n', mfilename) ...
    ,logFilePath,true);


%% Preallocation
if ~exist('nearestNeighbourCases','var')
    logMessage('Loading old results to continue with simulation.' ...
        ,logFilePath,false);
    if exist(configuration.oldResultsFileToLoad,'file')
        oldResults = load(configuration.oldResultsFileToLoad);
        additionalInformation = sprintf("These results are a " ...
            + "continuation of %s",configuration.oldResultsFileToLoad);
    else
        error("Old results file does not exist");
    end
    
    nearestNeighbourCases = getValuesFromStringEnumeration( ...
        configuration.nearestNeighbourCases,";","numeric");
    nearestNeighbourCases = sort(nearestNeighbourCases,'descend');
    maxNearestNeighbours = nearestNeighbourCases(1);
    
    theta = configuration.fibreAnglesTheta;
    yAxis = [0 1 0];
    
    phi = configuration.fibreAnglesPhi;
    zAxis = [0 0 1];
    
    logMessage("Comparing old with new configuration.",logFilePath,false);
    checkForEqualConfigurationsInNewAndOldSimulation( ...
        fileName,oldResults.configuration.fileName);
    
    checkForEqualConfigurationsInNewAndOldSimulation( ...
        nearestNeighbourCases ...
        ,sort(oldResults.nearestNeighbourCases,'descend'));
    
    checkForEqualConfigurationsInNewAndOldSimulation( ...
        averagingRegionForSpectralDensity ...
        ,oldResults.averagingRegionForSpecDens);
    
    checkForEqualConfigurationsInNewAndOldSimulation(theta ...
        ,oldResults.configuration.fibreAnglesTheta);
    checkForEqualConfigurationsInNewAndOldSimulation(phi ...
        ,oldResults.configuration.fibreAnglesPhi);
    
    checkForEqualConfigurationsInNewAndOldSimulation(omega0 ...
        ,oldResults.omega0);
    
    
    if ~(atomsToCalculate > oldResults.atomCounter)
        error("Atoms to calculate (%i) must be larger than the " ...
            + " atom counter`(%i).",atomsToCalculate ...
            ,oldResults.atomCounter);
    end
    logMessage("All good.",logFilePath,false);
    
    
    logMessage("Preallocation of arrays.",logFilePath,true);
    nearestNeighboursIDs = zeros(1,maxNearestNeighbours);
    nearestNeighbourDistancesPow3 = zeros(maxNearestNeighbours ...
        ,timeSteps,'single');
    
    polarAngle = zeros(maxNearestNeighbours,timeSteps,'single');
    azimuthAngle = zeros(maxNearestNeighbours,timeSteps,'single');
    
    sphericalHarmonicZerothOrder = zeros(maxNearestNeighbours,timeSteps ...
        ,'like',double(1j));
    sphericalHarmonicFirstOrder = zeros(maxNearestNeighbours,timeSteps ...
        ,'like',double(1j));
    sphericalHarmonicSecondOrder = zeros(maxNearestNeighbours,timeSteps ...
        ,'like',double(1j));
    
    corrFuncZerothOrder = zeros(1,timeSteps,'like',double(1j));
    corrFuncFirstOrder = zeros(1,timeSteps,'like',double(1j));
    corrFuncSecondOrder = zeros(1,timeSteps,'like',double(1j));
    
    logMessage("Preallocation ready.",logFilePath,true)
    
    logMessage("Set up variables based on old simulation." ...
        ,logFilePath,false);
    atomCounter = oldResults.atomCounter;
    r1Estimation = zeros(1,numberOfHs);
    r1Estimation(1:atomCounter) = oldResults.r1Estimation(1:atomCounter);
    calculatedAtomIndices = zeros(1,atomsToCalculate);
    calculatedAtomIndices(1:atomCounter) = ...
        oldResults.calculatedAtomIndices(1:atomCounter);
    randomSequenceOfAtoms = oldResults.randomSequenceOfAtoms;
    
    sumCorrFuncZerothOrder = oldResults.sumCorrFuncZerothOrder;
    sumCorrFuncFirstOrder = oldResults.sumCorrFuncFirstOrder;
    sumCorrFuncSecondOrder = oldResults.sumCorrFuncSecondOrder;
    
    corrFuncDirectories = createCorrFuncDirectoriesIfNotExist( ...
        resultsDirectory,gromacsSimulationDate);
    atomCounter = atomCounter + 1;
    logMessage("Created directories for correlation functions" ...
        + " if necessary.",logFilePath,false);
    
    logMessage(sprintf("This log file appends the log file under %s" ...
        ,oldResults.logFilePath),logFilePath,false);
    clear oldResults
    return;
end
meanPositions = single([mean(trajectoryX,2) mean(trajectoryY,2) ...
        mean(trajectoryZ,2)]);
atomTimer = zeros(1,atomsToCalculate);

%% Start simulation
logMessage(sprintf("You have chosen theta = %.4f and phi = %.4f" ...
    ,theta,phi),logFilePath,false);
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
        = calculateSpectralDensities( ...
        sumCorrFuncFirstOrder(1,:)/atomCounter ...
        ,sumCorrFuncSecondOrder(1,:)/atomCounter,omega0,deltaTInS ...
        ,averagingRegionForSpectralDensity);
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
            "   => Estimated relaxation rate for NN = %i " ...
            + "calculated atoms: %.4f",nearestNeighbourCases(1) ...
            ,r1Estimation(atomCounter)),logFilePath,false);
    end
    
    printEqualSignBreakLineToLogFile(logFilePath);
    atomCounter = atomCounter + 1;
end



