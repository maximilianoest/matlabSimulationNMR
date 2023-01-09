clc; clear all; close all; fclose('all');
% Description:
% In the orientation-dependency paper the authors said that in the
% cross-region is where both densities of solid lipid and lipid water H are
% larger than 0. They found that this region are 32% of the membran
% associated hydrogen nuceli. Here I am interested in the part of the lipid
% water pool that is interacting with the solid lipid part. Therefore, I
% divide the density of the interacting part by the non-interacting part.

%% add directories
addpath(genpath(sprintf('..%s..%slibrary',filesep,filesep)));
addpath(genpath(sprintf('..%s..%stxtFiles',createFilesepStringArray(2))));
constants = readConstantsFile('constants.txt');

%% check if paths are correct
filesDirectory = ...
   ['C:\Users\maxoe\Google Drive\Promotion\Simulation\RESULTS' ...
   '\densityDistributions\Plots\'];
fileNames = ["densityDistribution_DOPS" ...
    "densityDistribution_PLPC" "densityDistribution_PSM"];
dateToLoad = "20221122";
fileNames = dateToLoad + fileNames;
for file = fileNames
    filePath = sprintf("%s%s.mat",filesDirectory,file);
    if ~exist(filePath,'file')
        error('File cannot be found %s.',filePath);
    end
end
saving = 1;

%% load data

fig = initializeFigure();
fourthSubplotText = "";
for fileNr = 1:length(fileNames)
   lipidName = strrep(fileNames(fileNr) ...
       ,dateToLoad + "densityDistribution_","");
   data = load(sprintf("%s%s.mat",filesDirectory,fileNames(fileNr)));
   if lipidName ~= data.lipidName
       error("not the same lipide");
   end
   subFig(fileNr) = initializeSubplot(fig,2,2,fileNr); %#ok<SAGROW>
   title(sprintf("Density H-Atoms (Lipid: %s)", lipidName));
   legendEntries = {};
   avgLocations = data.avgLocations;
   plts = [];
   density = [];
   lowerIndices = [];
   upperIndices = [];
   moleculesToCompare = data.moleculesToCompare;
   for moleculeNr = 1:length(moleculesToCompare)
       densDistr = squeeze(mean(data.densityDistributionOnlyHdydrogen( ...
           :,moleculeNr,:),1));
       plts(end+1) = plot(avgLocations,densDistr); %#ok<SAGROW>
       if data.moleculesToCompare(moleculeNr) == "SOL"
           legendEntries{end+1} = "Water"; %#ok<SAGROW>
           coherentRegion = findCoherentRegion(find(densDistr > 0.1));
       else
           legendEntries{end+1} = "Lipid"; %#ok<SAGROW>
           coherentRegion = findCoherentRegion(find(densDistr < 0.1));
       end
       
       lowerBorder = avgLocations(coherentRegion(1));
       upperBorder = avgLocations(coherentRegion(end));
       
       density(moleculeNr,:) = densDistr'; %#ok<SAGROW>
       lowerIndices = sort([lowerIndices coherentRegion(1)]); 
       upperIndices = sort([upperIndices coherentRegion(end)]);
       
       plotVerticalLine(lowerBorder,120);
       plotVerticalLine(upperBorder,120);
   end
   
   hydrogenAtomWeight = data.atomWeightCollection( ...
       data.atomNameCollection == "H");
   averageDVolume = mean(data.dVolumeForTimeStep);
   
   fprintf("Lipid %s: \n",lipidName);
   %% water interaction fraction based on density distribution
   fprintf("   Water pool based on denstiy distribution: \n");
   numberOfHsInWaterPoolBasedOnDensity = ...
       density(moleculesToCompare == "SOL",:)/hydrogenAtomWeight ...
       *averageDVolume;
   fprintf("   Found %.2f hydrogen atoms in pool. \n" ...
       ,sum(numberOfHsInWaterPoolBasedOnDensity));
   numberOfHsInInteractionRegionOfWater = ...
       [numberOfHsInWaterPoolBasedOnDensity(1:lowerIndices(2)) ...
       numberOfHsInWaterPoolBasedOnDensity(upperIndices(1):end)];
   fprintf("   Found %.2f hydrogen atoms in interaction region. \n" ...
       ,sum(numberOfHsInInteractionRegionOfWater));
   numberOfHsInFreeWater = ...
       numberOfHsInWaterPoolBasedOnDensity(lowerIndices(2)+1:upperIndices(1)-1);
   fprintf("   Found %.2f hydrogen atoms in the free water region.\n" ...
       ,sum(numberOfHsInFreeWater));
   interactingWaterFraction = ...
       sum(numberOfHsInInteractionRegionOfWater) ...
       /sum(numberOfHsInWaterPoolBasedOnDensity);
   fprintf("   The water interaction fraction based on density is %.4f. \n" ...
       ,interactingWaterFraction)
   fprintf(" ------------------ \n");
   
   %% lipid interaction fraction based on density distribution
   fprintf("   Lipid pool based on density distribution: \n");
   numberOfHsInLipidPoolBasedOnDensity = ...
       density(moleculesToCompare == lipidName,:)/hydrogenAtomWeight ...
       *averageDVolume;
   fprintf("   Found %.2f hydrogen atoms in pool. \n" ...
       ,sum(numberOfHsInLipidPoolBasedOnDensity));
   numberOfHsInInteractionRegionOfLipid = ...
       numberOfHsInLipidPoolBasedOnDensity(lowerIndices(1):upperIndices(2));
   fprintf("   Found %.2f hydrogen atoms in interaction region. \n" ...
       ,sum(numberOfHsInInteractionRegionOfLipid));
   numberOfHsInFreeLipid = ...
       [numberOfHsInLipidPoolBasedOnDensity(1:lowerIndices(1)-1) ...
       numberOfHsInLipidPoolBasedOnDensity(upperIndices(2)+1:end)];
   fprintf("   Found %.2f hydrogen atoms in free lipid region.\n" ...
       ,sum(numberOfHsInFreeLipid)); 
   interactingLipidFraction = ...
       sum(numberOfHsInInteractionRegionOfLipid) ...
       /sum(numberOfHsInLipidPoolBasedOnDensity);
   fprintf("   The lipid interaction fraction based on density is %.4f. \n" ...
       ,interactingLipidFraction);
   fprintf(" ------------------ \n");
   poolFractionFactorEstimate = ...
       sum(numberOfHsInLipidPoolBasedOnDensity) ...
       / sum(numberOfHsInInteractionRegionOfWater);
   fprintf("   Estimate for pool fraction factor: %.4f \n" ...
       ,poolFractionFactorEstimate);
   
   fprintf(" ==============================\n");
   
   %% --------
   
   data.interactingWaterFraction = interactingWaterFraction;
   data.interactingLipidFraction = interactingLipidFraction;
   data.numberOfHsInLipidPool = numberOfHsInLipidPoolBasedOnDensity;
   data.numberOfHsInWaterPool = numberOfHsInWaterPoolBasedOnDensity;
   if saving
       save(sprintf("%s%s_withInteractionFractions.mat" ...
           ,filesDirectory,fileNames(fileNr)),'-struct','data');
   end
   fourthSubplotText = sprintf("%sLipid: %s\n  water interaction " ...
       +"fraction: %.4f \n  lipid interaction fraction: %.4f \n" ...
       +"  estimate for pool fraction factor: %.4f \n\n" ...
       ,fourthSubplotText,lipidName,interactingWaterFraction ...
       ,interactingLipidFraction,poolFractionFactorEstimate);
   
   ylabel("Density $[kg/m^3]$");
   xlabel("Location $[m]$");
   xticks(-5e-9 : 1e-9 : 5e-9);
   lgd = legend(plts,legendEntries);
   if fileNr == 3
       lgd.Location = "North";
   else
       lgd.Visible = "Off";
   end 
   
   
end
%%
ax = initializeSubplot(fig,2,2,fileNr+1);
text(0.05,0.5,fourthSubplotText,'interpreter','latex');
set(ax,'visible','off');
lgd = legend();
lgd.Visible = 'Off';
if saving
    saveFigureTo(filesDirectory ...
        ,"LocationDistributions_InteractionRegions","allLipids"...
        ,"onlyHydrogen");
end
    

function plt = plotVerticalLine(xPos,yValue)
plt = plot([xPos xPos],[0 yValue],'--k','LineWidth',0.7);
end

function indices = findCoherentRegion(indices)
differenceInIndices = diff(indices);
lengthOfCoherentRegion = 1;
for counter = 1:length(differenceInIndices)
    if differenceInIndices(counter) == 1
        lengthOfCoherentRegion(end) = lengthOfCoherentRegion(end) + 1;
    else
        lengthOfCoherentRegion(end+1) = 1; %#ok<AGROW>
    end
end
[maxLength,index] = max(lengthOfCoherentRegion);
if index == 1
    indices = indices(1:maxLength);
else
    startIndex = sum(lengthOfCoherentRegion(1:index-1)) + 1;
    endIndex = startIndex -1 + maxLength;
    indices = indices(startIndex : endIndex);
end
end

