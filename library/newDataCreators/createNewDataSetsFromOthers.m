function createNewDataSetsFromOthers()

%% load data
fileDirectories = [ ...
    "/daten/a/Relaxation/Lipids/simulation_max/DOPS/"...
    "/daten/a/Relaxation/Lipids/simulation_max/PLPC/" ...
    "/daten/a/Relaxation/Lipids/simulation_max/PSM/"];

fileNames = [ ...
    "20220110_DOPS_TIP4_Bilayer_50water_lipid_H_whole_dt2ps_simTime1000ns.mat" ...
    "20220401_PLPC_TIP4_Bilayer_50water_lipid_H_whole_dt2ps_simTime1000ns.mat" ...
    "20220509_PSM_TIP4_Bilayer_50water_lipid_H_whole_dt2ps_simTime1000ns.mat"];

filePaths = fileDirectories + fileNames;
timeStepSkips = [1 2 3 5 10];
shortings = [1];

%% save data as with other simulation parameters
for fileNr = 1:length(fileNames)
    fprintf("\n\n");
    filePath = filePaths(fileNr);
    fileName = fileNames(fileNr);
    fileDirectory = fileDirectories(fileNr);
    fprintf('Starting with %s \n',fileName);
    splittedFileName = strsplit(fileName,'_');
    gromacsSimulationDate = splittedFileName(1);
    whichLipid = splittedFileName(2);
    waterModel = splittedFileName(3);
    layerFormat = splittedFileName(4);
    waterCount = splittedFileName(5);
    constituent = splittedFileName(6);
    atom = splittedFileName(7);
    composingMode = splittedFileName(8);
    deltaT = replace(splittedFileName(9),{'dt','ps'},'');
    if strcmp(deltaT{1}(1),'0')
        afterComma = str2num(deltaT{1}(2:end));
        deltaT = afterComma/10;
    else
        deltaT = str2num(deltaT);
    end
    simulationDuration = str2num(replace(splittedFileName(10) ...
        ,{'simTime','ns.mat'},''));
    
    fprintf("deltaT: %dps, simulationDuration: %d ns \n",deltaT,simulationDuration);
    matObject = matfile(filePath);
    configuration = matObject.configuration;
    
    fprintf('Start loading data \n');
    load(filePath,'trajectories');
    wholeTrajectory = trajectories;
    clear trajectories
    [atomsCount,dimensionsCount,timeSteps] = size(wholeTrajectory );
    simulationSteps = round(timeSteps*shortings);
    simDurations = simulationDuration*shortings;
    for simulationTimeNr = 1:length(simulationSteps)
        for timeStepSkipNr = 1:length(timeStepSkips)
            timeStepSkip = timeStepSkips(timeStepSkipNr);
            newEnd = simulationSteps(simulationTimeNr);
            fprintf("New deltaT: %dps, simulationDuration: %d ns \n" ...
                ,deltaT*timeStepSkip,simDurations(simulationTimeNr));
            trajectories(:,1,:) = squeeze(single(wholeTrajectory(:,1 ...
                ,1:timeStepSkip:newEnd)));
            trajectories(:,2,:) = squeeze(single(wholeTrajectory(:,2 ...
                ,1:timeStepSkip:newEnd)));
            trajectories(:,3,:) = squeeze(single(wholeTrajectory(:,3 ...
                ,1:timeStepSkip:newEnd)));
            
            timeStep = deltaT*timeStepSkip;
            if timeStep < 1
                timeStep = strrep(num2str(timeStep),'.','');
            else 
                timeStep = num2str(timeStep);
            end
            
            savingName = sprintf( ...
                '%s_%s_%s_%s_%s_%s_%s_%s_dt%sps_simTime%ins.mat' ...
                ,gromacsSimulationDate,whichLipid ...
                ,waterModel,layerFormat,waterCount,constituent,atom ...
                ,composingMode,timeStep,simDurations(simulationTimeNr));
            

            pathToSave = sprintf('%s%s',fileDirectory,savingName);
            save(pathToSave,'trajectories','configuration','-v7.3');
            fprintf('Saved %s in directory %s \n',savingName,fileDirectory);
            clear trajectories
        end
    end
end

end