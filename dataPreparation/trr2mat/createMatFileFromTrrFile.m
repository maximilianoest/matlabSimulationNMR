clc
clear all

path2LipidData = "/daten/a/Relaxation/Lipids/simulation_max/"; 
directories = path2LipidData ...
    + ["DOPS/GROMACS/20220110_DOPS_TIP4_Monolayer_50water/" ...
    "PLPC/GROMACS/20220804_PLPC_TIP4_Monolayer_50water/" ...
    "PSM/GROMACS/20220804_PSM_TIP4_Monolayer_50water/"];

%% ATTENTION
% start and end time point given at seconds that should be loaded from the
% binary file that is saved from the .trr file.
startTimePoint = 40e-9;
endTimePoint = 50e-9;
 
%%
trrFileName = "prd";
path2BinFile = [pwd '/'];
simulationConfigurationFile = 'step7_production.mdp';
timeStepSkips = [100];
partialSimulationLengths = [1];



filePaths = string([]);
for directory = directories
    filePaths(end+1) = sprintf('%s%s.trr',directory,trrFileName); ...
        %#ok<SAGROW>
    if ~isfile(filePaths(end))
        error('%s does not exist',filePaths(end));
    end
end
fprintf('All files exist \n')

for filePath = filePaths
    splittedPath = strsplit(filePath,'/');
    lipid = splittedPath{end-3};
    simulationFileName = splittedPath{end-1};
    configFileId = fopen(sprintf('%s/%s',directory, ...
        simulationConfigurationFile),'r');
    configData = textscan(configFileId, '%s %s','Delimiter','=');
    
    configuration = {};
    for configElement = 1:size(configData{1},1)
        configurationName = configData{1}{configElement};
        configurationName = strrep(configurationName,' ','');
        configurationName = strrep(configurationName,'-','_');
        if configurationName ~= ';'
            try
                configurationInformation = configData{2}{configElement};
                configurationName = strrep(configurationName,'-','_');
                if isnan(str2double(configurationInformation))
                    configuration.(configurationName) = ...
                        configurationInformation;
                else
                    configuration.(configurationName) = ...
                        str2double(configurationInformation);
                end
                
            catch e
                warning(['For the configuration variable "' ...
                    configurationName '" no data are given.'
                    ' Take a look at the configuration file.']);
            end
        end
    end
    
    disp(['Lipid: ' lipid])
    disp(['Overall simulation steps: ' num2str(configuration.nsteps)])
    disp(['Gromacs simulation time step [ps]: ' num2str(configuration.dt)])
    disp(['Saving interval [#simulationSteps]: ' ...
        num2str(configuration.nstxout)])
    framesCount = configuration.nsteps/configuration.nstxout;
    timeBetweenDataPoints = configuration.dt * configuration.nstxout;
    disp(['Consequently, there are ' num2str(framesCount) ' frames ' ...
        'with a dT of ' num2str(timeBetweenDataPoints) ' ps.']);
    disp('--------------------------');
    
    startFrameSimulation = round( ...
        startTimePoint / (timeBetweenDataPoints*1e-12)) ;
    endFrameSimulation = round( ...
        endTimePoint / (timeBetweenDataPoints*1e-12));
    if endFrameSimulation > framesCount
        error("End frame was chosen too large.");
    end
    if startFrameSimulation > endFrameSimulation
        error("Start frame is larger than end frame.");
    end
    configuration.matFileSavedAtStartTimePointInS = startTimePoint;
    configuration.matFileSavedAtEndTimePointInS = endTimePoint;
    
    
    simulationTimeInPs = configuration.dt * configuration.nsteps;
    simulationTimeInPs_control = framesCount * timeBetweenDataPoints;
    if simulationTimeInPs ~= simulationTimeInPs_control
        error('simulation time is not equal');
    end
    
    disp(["Simulation originally was performed at :"]);
    disp(['Simulation duration [ps]: ' num2str(simulationTimeInPs)])
    disp(['Simulation duration [ns]: ' num2str(simulationTimeInPs/1000)])
    disp(['Simulation duration [us]: ' num2str(simulationTimeInPs/1000000)])
    
    disp(['Simulation will be read from binary as :']);
    disp(['Starting at time ' num2str(startTimePoint)  ...
        'sec and ending at time ' num2str(endTimePoint) ' sec']);
    cuttedSimulationTimeInPs = (endFrameSimulation - startFrameSimulation) ...
        * timeBetweenDataPoints;
    disp(['Thus, the available simulation duration is ' ...
        num2str(cuttedSimulationTimeInPs) ' ps']);
    disp('--------------------------');
    
    disp('The follwing simulation configurations will be saved ');
    disp(' sim. Dur. | deltaT :');
    
    reducedAndCuttedSimulationTimesInNs = ...
        cuttedSimulationTimeInPs * partialSimulationLengths/1000;
    
    for timeStepCountNr = 1:length(reducedAndCuttedSimulationTimesInNs)
        savingSimTime = strrep(num2str( ...
            reducedAndCuttedSimulationTimesInNs(timeStepCountNr)),'.','');
        for timeStepSkip = timeStepSkips       
            savingTimeBetweenDataPoints = configuration.dt  ...
                * configuration.nstxout * timeStepSkip;
            disp(['  ' num2str(savingSimTime) ' ns   |   ' ...
                num2str(savingTimeBetweenDataPoints) ' ps']);
        end
    end
    
    
    
    %% ATTENTION
    % x and z is changed so that the nerve fibre is oriented parallel to B0
    trr2matlab(filePath,'x',1250000);
    
    %% 
    tic
    [xdata] = readGmx2Matlab([path2BinFile 'xdata.binary'] ...
        ,startFrameSimulation,endFrameSimulation);
    fprintf('Time needed for processing binary file: %f s \n',toc);
    loadedFramesCount = size(xdata.trajectory,3);
    disp(['The number of loaded frames is: ' num2str(loadedFramesCount) ...
        ' which is equal to a simulation duration of ' ...
        num2str(loadedFramesCount*timeBetweenDataPoints) ' ns']);
    
    reducedLoadedFramesCount = round(loadedFramesCount ...
        * partialSimulationLengths);
    
    for timeStepCountNr = 1:length(reducedLoadedFramesCount)
        savingSimTime = strrep(num2str( ...
            reducedAndCuttedSimulationTimesInNs(timeStepCountNr)),'.','');
        for timeStepSkip = timeStepSkips
            trajectories = zeros(size(xdata.trajectory(:,:, ...
                1:timeStepSkip:reducedLoadedFramesCount(timeStepCountNr))) ...
                ,'single');
            trajectories(:,1,:) = single(xdata.trajectory(:,3, ...
                1:timeStepSkip:reducedLoadedFramesCount(timeStepCountNr)));
            trajectories(:,2,:) = single(xdata.trajectory(:,2, ...
                1:timeStepSkip:reducedLoadedFramesCount(timeStepCountNr)));
            trajectories(:,3,:) = single(xdata.trajectory(:,1, ...
                1:timeStepSkip:reducedLoadedFramesCount(timeStepCountNr)));
            
            savingTimeBetweenDataPoints = strrep(num2str(configuration.dt  ...
                * configuration.nstxout * timeStepSkip),'.','');
            savingName = sprintf('%s_%s_dt%sps_simTime%sns', ...
                simulationFileName,trrFileName ...
                ,savingTimeBetweenDataPoints,savingSimTime);
            
            tic
            pathToSave = sprintf('%s%s/%s',path2LipidData,lipid,savingName);
            save(pathToSave,'trajectories','configuration','-v7.3');
            fprintf('Saved %s in directory %s \n',trrFileName,pathToSave);
            fprintf('Time needed for saving .mat file: %f s \n', toc);
        end
    end
    
    clear trajectories xdata
    delete([path2BinFile '/xdata.binary']);
end









