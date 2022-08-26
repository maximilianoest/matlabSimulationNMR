function [trajectoryX,trajectoryY,trajectoryZ ...
    ,gromacsSimulationConfiguration] = loadTrajectoriesAndSimConfig( ...
    dataFilePath)

matObject = matfile(dataFilePath);
gromacsSimulationConfiguration = matObject.configuration;

load(dataFilePath,'trajectories');

trajectoryX = squeeze(single(trajectories(:,1,:))); %#ok<IDISVAR,NODEF>
trajectoryY = squeeze(single(trajectories(:,2,:)));
trajectoryZ = squeeze(single(trajectories(:,3,:)));

end