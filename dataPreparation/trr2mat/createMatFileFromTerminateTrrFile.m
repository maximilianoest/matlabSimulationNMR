clc; clear all; close all;

matFilesFolder = "/daten/a/Relaxation/MYELIN/Monolayer/GROMACS/";
gromacsFolderName =  ...
    "20230404_MYELIN_TIP4_50Water_ShortSimDur";
path2SimData = matFilesFolder + gromacsFolderName + "/";
    

trrFileName = "prd";
path2BinFile = [pwd '/'];
simulationConfigurationFile = 'step7_production.mdp';

filePath = path2SimData + trrFileName + ".trr";
if ~isfile(filePath)
    error("File does not exist.");
else
    fprintf("File exist.\n");
    
end

%% ATTENTION
% The simulation data loaded here are from a terminated simulation. Thus,
% the configuration data and the simulation data do not coincide.
% Therefore, information about the simulation length and the simulation
% deltaT are written directly into the trajectory file
 
%% 
configFileId = fopen(sprintf('%s%s',path2SimData ...
    ,simulationConfigurationFile),'r');
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

deltaTInPs = configuration.dt * configuration.nstxout;
deltaTInS = deltaTInPs * 1e-12;
fprintf("DeltaT = %.4d \n",deltaTInS);

%% ATTENTION
% x and z is changed so that the nerve fibre is oriented parallel to B0
tic
trr2matlab(filePath,'x',1250000);
fprintf("Time needed for creating binary file: %.5f.\n",toc);

%% create trajectories
tic
[xdata] = readGmx2Matlab(path2BinFile + "xdata.binary");
fprintf("Time needed for processing binary file: %fs \n",toc);

loadedFramesCount = size(xdata.trajectory,3);
configuration.nsteps = loadedFramesCount / 0.002;
fprintf("Number of loaded frames: %i",loadedFramesCount);
simDurInS = (loadedFramesCount - 1) * deltaTInS;
simDurInNs = simDurInS * 1e9;
fprintf("Simulation duration = %.4d\n",simDurInS);
numberOfAtoms = size(xdata.trajectory,1);
fprintf("Found %i atoms in the data set.\n",numberOfAtoms);

trajectories = zeros(numberOfAtoms,3,loadedFramesCount);
trajectories(:,1,:) = xdata.trajectory(:,3,:);
trajectories(:,2,:) = xdata.trajectory(:,2,:);
trajectories(:,3,:) = xdata.trajectory(:,1,:);

%% save data
deltaTForSavingName = strrep(num2str(deltaTInPs),'.','_');
simDurForSavingName = strrep(num2str(simDurInNs),'.','_');

savingName = sprintf("%s_%s_dt%sps_simDur%sns",gromacsFolderName ...
    ,trrFileName,deltaTForSavingName,simDurForSavingName);

pathToSave = sprintf('%s%s',matFilesFolder,savingName);

save(pathToSave,'trajectories','configuration','deltaTInS' ...
    ,'simDurInS','-v7.3');














