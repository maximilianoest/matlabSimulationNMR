%% script start
logMessage(sprintf('Starting external script: %s.m', mfilename),logFilePath ...
    ,true);

lags = round(configuration.fractionForLags*timeSteps);
logMessage(sprintf(['The lag is set to %d time steps, which ' ...
    'is equivalent to %.2f %%. This configuration only shortens the '...
    'correlation functions and NOT the simulation time.'],lags ...
    ,(configuration.fractionForLags)*100),logFilePath,false);



nearestNeighbours = configuration.nearestNeighbours;
if nearestNeighbours >= numberOfHs
    logMessage(['The number of nearest neighbours is higher than '
        'the number of possible atoms. PLEASE CHECK YOUR CONFIG ' ...
        'FILE!'],path2LogFile);
    error(['The number of nearest neighbours is higher than the ' ...
        'number of possible atoms. Please check your config file!']);
end

logMessage(sprintf(['Analysing %.f nearst neighbours of ' ...
    'overall %.f hydrogen atoms'],nearestNeighbours,numberOfHs) ...
    ,path2LogFile,false);



fibreAnglesTheta = deg2rad(getValuesFromStringEnumeration( ...
    configuration.fibreOrientations,';','numeric'));
fibreAnglesThetaCount = size(fibreAnglesTheta,2);
logMessage(['Found the orientations' sprintf(' %.f',rad2deg( ...
    fibreAnglesTheta))],path2LogFile,false);

fibreAnglesPhi = deg2rad(getValuesFromStringEnumeration( ...
    configuration.myelinPositions,';','numeric'));
fibreAnglesPhiCount = size(fibreAnglesPhi,2);
logMessage(['Found the positions' sprintf(' %.f',rad2deg( ...
    fibreAnglesPhi))],path2LogFile,false);