clear all; close all; clc; fclose('all');

addpath(genpath("..\..\library\"));

path2Data = "C:\Users\maxoe\Documents\Gromacs\testDatasets\MYELIN\Bilayer\";
fileName = "20221222_MYELIN_TIP4_Bilayer_50water_solidMyelin_H_whole_dt20ps_simTime200ns.mat";

fprintf("Loading data. \n");
load(path2Data + fileName);

%% plotting stuff
initializeFigure();
hold off
[~,~,timeSteps] = size(trajectories);
xlabel("X");
ylabel("Y");
zlabel("Z");

for timeStep = 1:timeSteps
   plot3(squeeze(trajectories(:,1,timeStep)) ...
       ,squeeze(trajectories(:,2,timeStep)) ...
       ,squeeze(trajectories(:,3,timeStep)),"*");
   pause(0.1);

end
