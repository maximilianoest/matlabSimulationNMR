clc
%% load data
fileDirectories = ["/daten/a/Relaxation/MYELIN/Monolaye/"];

fileNames = ["20230404_MYELIN_TIP4_50Water_ShortSimDur_prd_dt0_05ps_simDur0_7046ns"];

filePaths = fileDirectories + fileNames;
timeStepSkips = [40];
shortings = [1];

%% save data as with other simulation parameters
for fileNr = 1:length(fileNames)
    fprintf("\n\n");
    filePath = filePaths(fileNr);
    fileName = fileNames(fileNr);
    fileDirectory = fileDirectories(fileNr);
    disp(fileName);
    splittedFileName = strsplit(fileName,'_');
    gromacsSimulationDate = splittedFileName(1);
    whichLipid = splittedFileName(2);
    waterModel = splittedFileName(3);
    layerFormat = splittedFileName(4);
    waterCount = splittedFileName(5);
    constituent = splittedFileName(6);
    atom = "allAtoms";
    composingMode = splittedFileName(7);
    
    
    deltaTInPs = str2num(replace(splittedFileName(8),{'dt','ps'},''));
%     if strcmp(deltaT{1}(1),'0')
%         afterComma = str2num(deltaT{1}(2:end));
%         deltaT = afterComma/10;
%     else
%         deltaT = str2num(deltaT);
%     end
    simDurInNs = str2num(replace(splittedFileName(9) ...
        ,'simDur','') + "." +  replace(splittedFileName(10) ...
        ,'ns','')); %#ok<ST2NM>
    
    fprintf("deltaT: %dps, simulationDuration: %d ns \n",deltaTInPs ...
        ,simDurInNs);
    matObject = matfile(filePath);
    configuration = matObject.configuration;
    
    fprintf('Start loading data \n');
    load(filePath,'trajectories');
    wholeTrajectory = trajectories;
    clear trajectories
    [atomsCount,dimensionsCount,timeSteps] = size(wholeTrajectory );
    simulationTimes = round(timeSteps*shortings);
    simDurations = simDurInNs*shortings;
    for simulationTimeNr = 1:length(simulationTimes)
        for timeStepSkipNr = 1:length(timeStepSkips)
            timeStepSkip = timeStepSkips(timeStepSkipNr);
            newEnd = simulationTimes(simulationTimeNr);
            fprintf("deltaT: %dps, simulationDuration: %d ns \n" ...
                ,deltaTInPs*timeStepSkip,simDurations(simulationTimeNr));

            trajectories(:,1,:) = squeeze(single(wholeTrajectory(:,1,1:timeStepSkip:newEnd)));
            trajectories(:,2,:) = squeeze(single(wholeTrajectory(:,2,1:timeStepSkip:newEnd)));
            trajectories(:,3,:) = squeeze(single(wholeTrajectory(:,3,1:timeStepSkip:newEnd)));
            
            timeStep = deltaTInPs*timeStepSkip;
            if timeStep < 1
                timeStep = strrep(num2str(timeStep),'.','_');
            else 
                timeStep = num2str(timeStep);
            end
            
            simDurAsString = replace(num2str(simDurations( ...
                simulationTimeNr)),'.','_');
            
            savingName = sprintf( ...
                '%s_%s_%s_%s_%s_%s_%s_%s_dt%sps_simTime%sns.mat' ...
                ,gromacsSimulationDate,whichLipid ...
                ,waterModel,layerFormat,waterCount,constituent,atom ...
                ,composingMode,timeStep,simDurAsString);
            

            pathToSave = sprintf('%s%s',fileDirectory,savingName);
            save(pathToSave,'trajectories','configuration','-v7.3');
            fprintf('Saved %s in directory %s \n',savingName,fileDirectory);
            clear trajectories
        end
    end
end
