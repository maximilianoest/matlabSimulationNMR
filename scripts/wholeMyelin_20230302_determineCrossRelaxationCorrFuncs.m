logMessage(sprintf('Starting the external script: %s.m \n\n', mfilename) ...
    ,logFilePath,true);

if ~exist('nearestNeighbourCases','var')
    
    nearestNeighbourCases = getValuesFromStringEnumeration( ...
        configuration.nearestNeighbourCases,";","numeric");
    nearestNeighbourCases = sort(nearestNeighbourCases,'descend');
    maxNearestNeighbours = nearestNeighbourCases(1);
    atomCounter = 1;
    
    theta = configuration.fibreAnglesTheta;
    yAxis = [0 1 0];
    
    phi = configuration.fibreAnglesPhi;
    zAxis = [0 0 1];
    
    logMessage("Preallocation of arrays.",logFilePath,true);
    nearestNeighboursIDs = zeros(1,maxNearestNeighbours);
    nearestNeighbourDistancesPow3 = zeros(maxNearestNeighbours ...
        ,timeSteps,'single');
    
    polarAngle = zeros(maxNearestNeighbours,timeSteps,'single');
    azimuthAngle = zeros(maxNearestNeighbours,timeSteps,'single');
    
%     sphericalHarmonicZerothOrder = zeros(maxNearestNeighbours,timeSteps ...
%         ,'like',double(1j));
%     sphericalHarmonicFirstOrder = zeros(maxNearestNeighbours,timeSteps ...
%         ,'like',double(1j));
%     sphericalHarmonicSecondOrder = zeros(maxNearestNeighbours,timeSteps ...
%         ,'like',double(1j));
    sphericalHarmonic = zeros(3,maxNearestNeighbours,timeSteps ...
        ,'like',double(1j));
    
%     corrFuncZerothOrder = zeros(1,timeSteps,'like',double(1j));
%     corrFuncFirstOrder = zeros(1,timeSteps,'like',double(1j));
%     corrFuncSecondOrder = zeros(1,timeSteps,'like',double(1j));

    corrFunc = zeros(3,length(nearestNeighbourCases) ...
        ,timeSteps,'like',double(1j));
    
    sumCorrFuncZerothOrder = zeros(1,timeSteps,'like',double(1j));
    sumCorrFuncFirstOrder = zeros(1,timeSteps,'like',double(1j));
    sumCorrFuncSecondOrder = zeros(1,timeSteps,'like',double(1j));
    
    logMessage("Preallocation ready.",logFilePath,true)
    
    if ~configuration.reloadOldSimulation
        logMessage("No old simulation results will be loaded" ...
            ,logFilePath,false);
    elseif exist(configuration.oldResultsFileToLoad,'file')
        logMessage("Load old simulation results",logFilePath,false);
        oldResults = load(configuration.oldResultsFileToLoad);
        additionalInformation = sprintf("These results are a " ...
            + "continuation of %s",configuration.oldResultsFileToLoad);
        logMessage(additionalInformation,logFilePath,false);
        
        logMessage("Comparing old with new configuration values." ...
            ,logFilePath,false);
        checkForEqualConfigurationsInNewAndOldSimulation( ...
            fileName,oldResults.configuration.fileName);
        
        checkForEqualConfigurationsInNewAndOldSimulation( ...
            nearestNeighbourCases ...
            ,sort(oldResults.nearestNeighbourCases,'descend'));
        
        checkForEqualConfigurationsInNewAndOldSimulation( ...
            averagingRegionForSpectralDensity ...
            ,oldResults.averagingRegionForSpectralDensity);
        
        checkForEqualConfigurationsInNewAndOldSimulation(theta ...
            ,oldResults.configuration.fibreAnglesTheta);
        checkForEqualConfigurationsInNewAndOldSimulation(phi ...
            ,oldResults.configuration.fibreAnglesPhi);
        
        checkForEqualConfigurationsInNewAndOldSimulation(omega0 ...
            ,oldResults.omega0);
        
        logMessage("All good.",logFilePath,false);
        logMessage("Set up variables based on old simulation." ...
            ,logFilePath,false);
        atomCounter = oldResults.atomCounter;
        r1Estimation = zeros(1,numberOfHs);
        r1Auto = oldResults.r1Auto(1:atomCounter);
        r1Cross = oldResults.r1Cross(1:atomCounter);
        calculatedAtomIndices = zeros(1,atomsToCalculate);
        calculatedAtomIndices(1:atomCounter) = ...
            oldResults.calculatedAtomIndices(1:atomCounter);
        randomSequenceOfAtoms = oldResults.randomSequenceOfAtoms;
        atomTimer = zeros(1,atomsToCalculate);
        atomTimer(1:atomCounter) = oldResults.atomTimer(1:atomCounter);
        
        sumCorrFuncZerothOrder = oldResults.sumCorrFuncZerothOrder;
        sumCorrFuncFirstOrder = oldResults.sumCorrFuncFirstOrder;
        sumCorrFuncSecondOrder = oldResults.sumCorrFuncSecondOrder;
        
        atomCounter = atomCounter + 1;
        
        logMessage(sprintf("This log file appends the log file under %s" ...
            ,oldResults.logFilePath),logFilePath,false);
        clear oldResults
    else
        error("Old results file %s does not exist" ...
            ,configuration.oldResultsFileToLoad);
    end
    
    corrFuncDirectories = createCorrFuncDirectoriesIfNotExist( ...
        resultsDirectory,gromacsSimulationDate);
    logMessage("Created directories for correlation functions" ...
        + " if necessary.",logFilePath,false);
    
    groupsToSearchIn = getValuesFromStringEnumeration( ...
        configuration.groupsToSearchIn,";","string");
    membraneIndicesLine = find(groupsToSearchIn == "MEMB");
    waterIndicesLine = find(groupsToSearchIn == "Water");
    unassignedGroup = getValuesFromStringEnumeration( ...
        configuration.unassignedGroup,";","string");
    indexFilePath = sprintf("%s%s",dataDirectory ...
        ,configuration.indexFileName);
    matFileIndices = getMatFileIndicesFromIndexFile(indexFilePath ...
        ,groupsToSearchIn,unassignedGroup);
    
    randomSequenceOfMemebraneAtoms = matFileIndices{membraneIndicesLine} ...
        (randperm(length(matFileIndices{membraneIndicesLine})));
    
    return;
    
end

%% Start simulation
meanPositions = single([mean(trajectoryX,2) mean(trajectoryY,2) ...
        mean(trajectoryZ,2)]);

drawDistributionAccordingToMatFileIndices(meanPositions ...
    ,groupsToSearchIn,matFileIndices);
saveFigureTo(resultsDirectory,whichLipid,matlabSimulationDate ...
    ,"hydrogenAtomsDistributionAccordingToGroup",true);
close(gcf);

logMessage(sprintf("You have chosen theta = %.4f and phi = %.4f" ...
    ,theta,phi),logFilePath,false);
logMessage("Rotating dataset.",logFilePath);
rotationMatrixPhi = get3DRotationMatrix(deg2rad(phi),zAxis);
rotationMatrixTheta = get3DRotationMatrix(deg2rad(theta),yAxis);
totalRotationMatrix = rotationMatrixTheta*rotationMatrixPhi;
[trajectoryX,trajectoryY,trajectoryZ]  ...
    = rotateTrajectoriesWithRotationMatrix(totalRotationMatrix ...
    ,trajectoryX,trajectoryY,trajectoryZ);
for atomNumber = randomSequenceOfMemebraneAtoms( ...
        atomCounter:atomsToCalculate)
    
     atomTimerStart = tic;
    calculatedAtomIndices(atomCounter) = atomNumber; %#ok<SAGROW>
    
    [trajectoryX,trajectoryY,trajectoryZ] = calculateRelativePositions( ...
        trajectoryX,trajectoryY,trajectoryZ,atomNumber);
    [nearestNeighboursIDsInWaterPool,nearestNeighbourDistancesPow3] =  ...
        findNearestNeighboursIDs(maxNearestNeighbours ...
        ,trajectoryX(matFileIndices{waterIndicesLine},:) ...
        ,trajectoryY(matFileIndices{waterIndicesLine},:) ...
        ,trajectoryZ(matFileIndices{waterIndicesLine},:));
    nearestNeighboursIDs = matFileIndices{waterIndicesLine} ...
        (nearestNeighboursIDsInWaterPool);
    
    [polarAngle,azimuthAngle] = ...
        transformToSphericalCoordinates( ...
        trajectoryX(nearestNeighboursIDs,:) ...
        ,trajectoryY(nearestNeighboursIDs,:) ...
        ,trajectoryZ(nearestNeighboursIDs,:));            
    [sphericalHarmonic(1,:,:),sphericalHarmonic(2,:,:) ...
        ,sphericalHarmonic(3,:,:)] = calculateSphericalHarmonics( ...
        polarAngle,azimuthAngle,nearestNeighbourDistancesPow3);
    
%     corrFuncZerothOrder = ...
%         calculateCorrFuncForMultipleNNCasesAsMatrix( ...
%         sphericalHarmonicZerothOrder,nearestNeighbourCases);
%     corrFuncFirstOrder = ...
%         calculateCorrFuncForMultipleNNCasesAsMatrix( ...
%         sphericalHarmonicFirstOrder,nearestNeighbourCases);
%     corrFuncSecondOrder = ...
%         calculateCorrFuncForMultipleNNCasesAsMatrix( ...
%         sphericalHarmonicSecondOrder,nearestNeighbourCases);
    
    parfor i = 1:3
        corrFunc(i,:,:) = calculateCorrFuncForMultipleNNCasesAsMatrix( ...
            squeeze(sphericalHarmonic(i,:,:)),nearestNeighbourCases);
    end
    
    whichDirectory = resultsDirectory + corrFuncDirectories( ...
        corrFuncDirectories == gromacsSimulationDate + "_" ...
        + getNamesOfVariblesAsArray(corrFunc) + "ZerothOrder" + filesep);
    if ~isfolder(whichDirectory)
        error("Chosen wrong directory: %s",whichDirectory);
    end
    saveCorrFuncTo(whichDirectory,squeeze(corrFunc(1,:,:)),atomNumber ...
        ,nearestNeighbourCases,simulationDurationInNs,deltaTInPs ...
        ,theta,phi);
    
    whichDirectory = resultsDirectory + corrFuncDirectories( ...
        corrFuncDirectories == gromacsSimulationDate + "_" ...
        + getNamesOfVariblesAsArray(corrFunc) + "FirstOrder" + filesep);
    if ~isfolder(whichDirectory)
        error("Chosen wrong directory: %s",whichDirectory);
    end
    saveCorrFuncTo(whichDirectory,squeeze(corrFunc(2,:,:)),atomNumber ...
        ,nearestNeighbourCases,simulationDurationInNs,deltaTInPs...
        ,theta,phi);
    
    whichDirectory = resultsDirectory + corrFuncDirectories( ...
        corrFuncDirectories == gromacsSimulationDate + "_" ...
        + getNamesOfVariblesAsArray(corrFunc) + "SecondOrder"  + filesep);
    if ~isfolder(whichDirectory)
        error("Chosen wrong directory: %s",whichDirectory);
    end
    saveCorrFuncTo(whichDirectory,squeeze(corrFunc(3,:,:)),atomNumber ...
        ,nearestNeighbourCases,simulationDurationInNs,deltaTInPs...
        ,theta,phi);
    
    sumCorrFuncZerothOrder = sumCorrFuncZerothOrder ...
        + squeeze(corrFunc(1,:,:));
    sumCorrFuncFirstOrder = sumCorrFuncFirstOrder ...
        + squeeze(corrFunc(2,:,:));
    sumCorrFuncSecondOrder = sumCorrFuncSecondOrder ...
        + squeeze(corrFunc(3,:,:));
    
    [avgSpecDensZerothOrder,avgSpecDensFirstOrder ...
        ,avgSpecDensSecondOrder] ...
        = calculateSpecDensForZerothFirstAndSecondOrder( ...
        sumCorrFuncZerothOrder(1,:)/atomCounter ...
        ,sumCorrFuncFirstOrder(1,:)/atomCounter ...
        ,sumCorrFuncSecondOrder(1,:)/atomCounter, omega0,deltaTInS ...
        ,averagingRegionForSpectralDensity);
   
    r1Auto(atomCounter) = calculateR1AutoWithSpecDens( ...
        avgSpecDensZerothOrder,avgSpecDensFirstOrder ...
        ,avgSpecDensSecondOrder,dipolDipolConstant); %#ok<SAGROW>
    r1Cross(atomCounter) = calculateR1CrossWithSpecDens( ...
        avgSpecDensZerothOrder,avgSpecDensSecondOrder ...
        ,dipolDipolConstant); %#ok<SAGROW>
    
    logMessage(sprintf('Calculated %i atom(s)',atomCounter),logFilePath);
    atomTimer(atomCounter) = toc(atomTimerStart); %#ok<SAGROW>
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
        updateSimulationProgressForCrossAndAutoCorrelation( ...
            sumCorrFuncZerothOrder(1,:)/atomCounter ...
            ,sumCorrFuncFirstOrder(1,:)/atomCounter ...
            ,sumCorrFuncSecondOrder(1,:)/atomCounter ...
            ,avgSpecDensZerothOrder,avgSpecDensFirstOrder ...
            ,avgSpecDensSecondOrder,deltaTInS ...
            ,averagingRegionForSpectralDensity,r1Cross,r1Auto)
        
        
        saveFigureTo(resultsDirectory,whichLipid,matlabSimulationDate ...
            ,sprintf("simulationProgress_simDur%sns_dT%sps" ...
            ,getDecimalNumberAsSeparatedWithUnderScore( ...
            simulationDurationInNs) ...
            ,getDecimalNumberAsSeparatedWithUnderScore(deltaTInPs)),false);
        close(gcf);
        
        logMessage(sprintf( ...
            "   => Estimated auto-relaxation rate for NN = %i %.4f" ...
            ,nearestNeighbourCases(1),r1Auto(atomCounter)) ...
            ,logFilePath,false);
        logMessage(sprintf( ...
            "   => Estimated cross-relaxation rate for NN = %i %.4f" ...
            ,nearestNeighbourCases(1),r1Cross(atomCounter)) ...
            ,logFilePath,false);
    end
    
    printEqualSignBreakLineToLogFile(logFilePath);
    atomCounter = atomCounter + 1;
end





















