clc
clear all

path2Lipid = '/home/moesterreicher/GROMACS/DOPS/'; 
simulationDataFolder = '20220509_DOPS_TIP4_Bilayer_50water';
path2BinFile = [pwd '/'];
path2File = [path2Lipid simulationDataFolder '/'];
namesOfFiles = ["lipid_H_whole"];

for fileNr = 1:length(namesOfFiles)
    filePath = path2File + namesOfFiles(fileNr) + '.trr';
    if ~isfile(filePath)
       error('%s does not exist',filePath); 
    end
end
fprintf('All files exist \n')


simulationConfigurationFile = 'step7_production.mdp';
pathToSave = '';

lipid = strsplit(path2Lipid,'/');
lipid = lipid{end-1};

configFileId = fopen([path2Lipid '/' simulationDataFolder '/' ...
    simulationConfigurationFile],'r');
data = textscan(configFileId, '%s %s','Delimiter','=');

configuration = {};
for element = 1:size(data{1},1)
    configurationName = data{1}{element};
    configurationName = strrep(configurationName,' ','');
    configurationName = strrep(configurationName,'-','_');
    if configurationName ~= ';'
        try
            configurationInformation = data{2}{element};
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
disp(['Simulation time step [ps]: ' num2str(configuration.dt)])
disp(['Saving interval [#simulationSteps]: ' ...
    num2str(configuration.nstxout)])

timeBetweenDataPoints = configuration.dt * configuration.nstxout;
disp(['Time between two data points [ps]: ' ...
    num2str(timeBetweenDataPoints)])
timeBetweenDataPoints = strrep(num2str(timeBetweenDataPoints),'.','');

simulationTime = configuration.dt * configuration.nsteps;
disp(['Simulation time [ps]: ' num2str(simulationTime)])
disp(['Simulation time [ns]: ' num2str(simulationTime/1000)])
disp(['Simulation time [us]: ' num2str(simulationTime/1000000)])

simTime = strrep(num2str(simulationTime/1000),'.','');



for fileNameNr = 1: length(namesOfFiles)
    nameOfFile = namesOfFiles(fileNameNr);
    fileName = sprintf('%s_%s_dt%sps_simTime%sns', ...
        simulationDataFolder,nameOfFile,timeBetweenDataPoints,simTime);
    %% ATTENTION
    % x and z is changed so that the fibre is oriented parallel to B0
    trr2matlab(path2File + nameOfFile + '.trr','x',1);
    
    tic
    [xdata] = readGmx2Matlab([path2BinFile 'xdata.binary']);
    fprintf('Time needed for processing binary file: %f s \n',toc);
    
    trajectories=zeros(size(xdata.trajectory),'single');
    trajectories(:,1,:)=single(xdata.trajectory(:,3,:));
    trajectories(:,2,:)=single(xdata.trajectory(:,2,:));
    trajectories(:,3,:)=single(xdata.trajectory(:,1,:));
    
    delete([path2BinFile '/xdata.binary']);
    
    tic
    save([pathToSave fileName],'trajectories','configuration', '-v7.3');
    disp(['Saved ' fileName ' into current directory'])
    fprintf('Time needed for saving .mat file: %f s \n', toc);
    clear trajectories xdata
    
end








