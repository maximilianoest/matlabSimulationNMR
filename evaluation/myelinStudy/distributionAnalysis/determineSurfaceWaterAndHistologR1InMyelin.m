%%
clc; close all; clear all; fclose('all');


addpath(genpath(sprintf('..%s..%s..%slibrary',filesep,filesep,filesep)));
addpath(genpath(sprintf('..%s..%s..%stxtFiles' ...
    ,createFilesepStringArray(3))));
constants = readConstantsFile('constants.txt');

resultsDir = sprintf("..%s..%s..%sRESULTS%swholeMyelin" ...
    + "_densityDistributions%s",createFilesepStringArray(5));
if ~exist(resultsDir,'dir')
    warning('The results directory does not exist but will be created.');
    mkdir(resultsDir)
end

%% set up
windowSize = 11;
surfWaterDensBorder = 0.95;
borderShiftOffset = 0;
saving = 0;


%% density data
densDistribFileName = "20230405_myelinDistributionData_averaged";
densData = load(resultsDir + densDistribFileName + ".mat");

locations = densData.avgDiscreteLocations(:,2:end-1);
dVol = mean(densData.dVolumeForTimeStep);

waterIndex = (densData.moleculeNameCollection == "WATER");
membraneIndices = (densData.moleculeNameCollection ~= "WATER") ...
    .* (densData.moleculeNameCollection ~= "POT") ...
    .* (densData.moleculeNameCollection ~= "CLA");
hydrogenIndex = densData.atomNameCollection == "H";

densDistribWater = squeeze(sum( ...
    densData.avgDensity_locMolAtom(:,waterIndex,:),2));
densDistribMembr = squeeze(sum( ...
    densData.avgDensity_locMolAtom(:,~waterIndex,:),2));

modelWaterMass = sum(sum(densDistribWater)) * dVol;
modelMembrMass = sum(sum(densDistribMembr)) * dVol;

densDistribWaterH = meanFilter(sum( ...
    densData.avgDensity_locMolAtom(2:end-1,waterIndex,hydrogenIndex) ...
    ,2),windowSize);
densDistribMembrH = meanFilter(sum( ...
    densData.avgDensity_locMolAtom(2:end-1,~waterIndex,hydrogenIndex) ...
    ,2),windowSize);
modelWaterHMass = sum(sum(densDistribWaterH))*dVol;
modelMembrHMass = sum(sum(densDistribMembrH))*dVol;


%% Schyboll2019 water-mass ratio
% According to Schyboll2019, the water mass ratio is defined as:
% waterMassRatio = waterWeight/(overallMembrWeight+waterWeight)
% CHL, DPPCm POPE, GalC, GalS
% cWeight = (27*27+23*40+23*39+23*40+4*40) ...
%     * densData.atomWeightCollection(densData.atomNameCollection == "C");
% hWeight = (27*46+23*80+23*76+23*77+4*76) ...
%     * densData.atomWeightCollection(densData.atomNameCollection == "H");
% oWeight = (27*1+23*8+23*8+23*7+4*11) ...
%     * densData.atomWeightCollection(densData.atomNameCollection == "O");
% nWeight = (27*0+23*1+23*1+23*1+4*1) ...
%     * densData.atomWeightCollection(densData.atomNameCollection == "N");
% pWeight = (27*0+23*1+23*1+23*0+4*0) ...
%     * densData.atomWeightCollection(densData.atomNameCollection == "P");
% sWeight = (27*0+23*0+23*0+23*0+4*1) ...
%     * densData.atomWeightCollection(densData.atomNameCollection == "S");
% 
% overallMembrWeight = 2*(cWeight + hWeight + oWeight + nWeight + pWeight ...
%     + sWeight);
% 
% waterWeight = 4699*(2* ...
%     densData.atomWeightCollection(densData.atomNameCollection == "H") ...
%     + densData.atomWeightCollection(densData.atomNameCollection == "O"));
% 
% membraneHWeight = 2*hWeight;
% waterHWeight = 4699*2* ...
%     densData.atomWeightCollection(densData.atomNameCollection == "H");
% hWeightRatio = waterHWeight/(membraneHWeight+waterHWeight);
% 
% fprintf("<strong> In Schyboll2019 model: </strong> \n");
% waterMassRatioSchyboll2019 = waterWeight/(overallMembrWeight+waterWeight);
% fprintf("water mass ratio: %.4f\n",waterMassRatioSchyboll2019);
% fprintf("hydrogen mass ratio: %.4f\n",hWeightRatio);


%%
fprintf("<strong> In this model: </strong> \n");

waterMassRatio = modelWaterMass/(modelWaterMass+modelMembrMass);
fprintf("  Water mass: %.4d kg.\n  Membrane mass: %.4d kg.\n" ...
    ,modelWaterMass,modelMembrMass);
fprintf("  => water mass ratio: %.4f\n",waterMassRatio);

waterMassRatioH = modelWaterHMass/(modelWaterHMass+modelMembrHMass);
fprintf("  Water H mass: %.4d kg.\n  Membrane H mass: %.4d kg.\n" ...
    ,modelWaterHMass,modelMembrHMass);
fprintf("  => water mass ratio: %.4f\n",waterMassRatioH);

surfWaterBorderIndices = densDistribWaterH < ...
    surfWaterDensBorder * mean(densDistribWaterH(end/2-100:end/2+100));
surfWaterBorders = find(diff(surfWaterBorderIndices) ~= 0) ...
    + [borderShiftOffset, -borderShiftOffset];
surfWaterBorderIndices(1:surfWaterBorders(1)) = 1;
surfWaterBorderIndices(surfWaterBorders(2):end) = 1;
if length(surfWaterBorders) ~= 2
    error("The borders for surface water are not chosen correctly.");
end

densSurfWaterH = nan(1,length(densDistribWaterH));
densSurfWaterH(surfWaterBorderIndices) = densDistribWaterH( ...
    surfWaterBorderIndices);

initializeFigure();
plot(locations,densDistribWaterH);
plot(locations,densSurfWaterH);
plot(locations,densDistribMembrH);
plotVerticalLine(locations(surfWaterBorders(1)),max(densDistribWaterH));
plotVerticalLine(locations(surfWaterBorders(2)),max(densDistribWaterH));
legend("Free H$_{2}$O","Surf. H$_{2}$O","Membrane",'Location','south');
title("Hydrogen atoms density");
xlabel('x position $[m]$');
ylabel('density $[kg/m^3]$');

%% relaxation stuff

r1DataFolder = "C:\Users\maxoe\Google Drive\Promotion\Simulation" ...
    + "\RESULTS\wholeMyelin_relaxationRates\";
r1DataFileName = "20230419_solidMyelinAndMyelinWater_CompartmentAndCrossR1";
r1Data = load(r1DataFolder + r1DataFileName + ".mat");
fieldStrength3TeslaIndex = r1Data.fieldStrengths < 3.00001 ...
    & r1Data.fieldStrengths > 2.99999;
fieldStrength1comma5TeslaIndex = r1Data.fieldStrengths < 1.50001 ...
    & r1Data.fieldStrengths > 1.499999;
fieldStrength7TeslaIndex = r1Data.fieldStrengths < 7.00001 ...
    & r1Data.fieldStrengths > 6.99999;

%% determine surface water R1

% based on Tsukiashi2018
% T1 vs. temperature: y = 0.058071*x + 1.7623 with x = Temperature in °C
temperature = 36.6; % °C from Obermayer
r1Data.r1_free =  1/(0.058071 * temperature + 1.7623);

fprintf("\n  Free water R1: %.4f \n", r1Data.r1_free);

% equation to determine lipid water R1 based on pool fraction factor:
% r1Eff_MW = f_surf * r1Surf_MW + f_free * r1_free
% ==> R1_surf = R1_tot/f_surf - R1_free*f_free/f_surf

frequencyWaterH = densData.avgFrequency_locMolAtom(:,waterIndex ...
    ,densData.atomNameCollection == "H");
overallWaterHCount = sum(frequencyWaterH);
surfaceWaterHCount = sum(frequencyWaterH(1:surfWaterBorders(1))) ...
    + sum(frequencyWaterH(surfWaterBorders(2):end));
f_surf_MD = surfaceWaterHCount / overallWaterHCount;
f_free_MD = 1-f_surf_MD;
fprintf("  f_surf = %.4f , f_free = %.4f\n",f_surf_MD,f_free_MD);

r1Data.r1Surf_MW = r1Data.r1Eff_MW/f_surf_MD ...
    - r1Data.r1_free*f_free_MD/f_surf_MD;
fprintf("  Surface water R1 at 3T: %.4f \n" ...
    ,r1Data.r1Surf_MW(fieldStrength3TeslaIndex));

waterMoleculeWeight = 2 * densData.atomWeightCollection( ...
    densData.atomNameCollection == "H") + densData.atomWeightCollection( ...
    densData.atomNameCollection == "O");
surfWaterMoleculeCount = surfaceWaterHCount/2;
surfWaterMass = surfWaterMoleculeCount * waterMoleculeWeight;
fprintf("\n  There are %.0f water molecules at the surface.\n" ...
    ,surfWaterMoleculeCount);
fprintf("  The surface water weight is: %.4d\n",surfWaterMass);

%% According to Morell1999 https://www.ncbi.nlm.nih.gov/books/NBK28221/
fprintf("\n<strong> Histological values based on literature: </strong> \n");
% R1 calculated with masses or with hydrogen density leads to the same
% results. Because the hydrogen density is essential for the signal, the
% calculations here are based on this density.
% water content is 40% -> 40% water, 60% non-water

waterContentInMyelin = 0.4; %according to Morell1999

%% According to studies in brainConstitution.m

constitutionData = load("C:\Users\maxoe\Google Drive\Promotion" ...
    + "\Simulation\RESULTS\wholeMyelin_brainConstitution" ...
    + "\20230419_wmAndGMCompositionBasedOnLiterature.mat");

wmWaterContent = constitutionData.wmWaterContent;
myelinWaterContent = constitutionData.wmWaterContent ...
    * constitutionData.myelinWaterContent;

% waterContentInMyelin = myelinWaterFraction/(myelinWaterFraction +
%                        myelinSolidFraction)
% -> solved for myelinSolidFraction: contains mass lipids and proteins
solidMyelinFraction = myelinWaterContent/waterContentInMyelin ...
    - myelinWaterContent;
fprintf("  Solid myelin fraction: %.4f\n",solidMyelinFraction);

% Thus, 0.0802 myelin water fraction comes hand in hand with 0.1203 solid
% myelin fraction. ESIVR is then the same like 0.4/0.6
effectiveSurfInteractionVolRatio = ...
    myelinWaterContent / solidMyelinFraction;
fprintf("  ESIVR for myelin: %.4f\n" ...
    ,effectiveSurfInteractionVolRatio);

% -> for given solid myelin mass, known from simulation and density data,
% the histological myelin water mass is the ESIVR multiplied by the
% given solid myelin mass of the lipids.
histMyelinWaterMass = modelMembrMass * effectiveSurfInteractionVolRatio;
histMyelinWaterHCount = histMyelinWaterMass / waterMoleculeWeight * 2;
fprintf("  histological myelin water H count based on lipid weight" ...
    + " in MD sim.: %.4f\n",histMyelinWaterHCount);

histFreeWaterHCount = histMyelinWaterHCount - surfaceWaterHCount;
fprintf("  free water H count in MD sim.: %.4f\n",histFreeWaterHCount);

histF_free = histFreeWaterHCount/histMyelinWaterHCount;
histF_surf = surfaceWaterHCount/histMyelinWaterHCount;
fprintf("  histolgical f_surf = %.4f, f_free = %.4f\n"  ...
    ,histF_surf,histF_free);

r1Data.r1Hist_MW = histF_free*r1Data.r1_free + histF_surf*r1Data.r1Surf_MW;
fprintf("  Healthy MW R1 at 1.5T: %.4f\n" ...
    ,r1Data.r1Hist_MW(fieldStrength1comma5TeslaIndex));
fprintf("  Healthy MW R1 at 3T: %.4f\n" ...
    ,r1Data.r1Hist_MW(fieldStrength3TeslaIndex));
fprintf("  Healthy MW R1 at 7T: %.4f\n" ...
    ,r1Data.r1Hist_MW(fieldStrength7TeslaIndex));

%% plotting
initializeFigure();

plot(r1Data.fieldStrengths,r1Data.r1Surf_MW);
plot(r1Data.fieldStrengths,r1Data.r1Hist_MW);
legend("R$_{1,surf}$", "R1$_{MW,hist}$");
xlabel("Field strength [Tesla]");
ylabel("R$_{1,MW}$ [Hz]");

if saving
    saveFigureTo(r1DataFolder,datestr(now,'yyyymmdd'),"myelinWater" ...
        ,"histologicalMyelinWaterR1",true);
end

fig = initializeFigure();
plot(r1Data.fieldStrengths,r1Data.r1Eff_SM);
plot(r1Data.fieldStrengths,r1Data.r1Hist_MW);

axis([0 inf 0 10]);
legend("R$^{eff}_{1,SM}$","R$^{eff}_{1,MW}$");

title("Histological compartment R$_1$");
xlabel("Field strength [Tesla]");
ylabel("R$_1^{eff}$");
xticks(0:0.5:r1Data.fieldStrengths(end));

if saving
    saveFigureTo(r1DataFolder,datestr(now,'yyyymmdd') ...
        ,"solidMyelinAndMyelinWater" ...
        ,"histologicalEffectiveCompartmentR1s",true);
    
end


if saving
    save(r1DataFolder + datestr(now,"yyyymmdd") ...
        +"_solidMyelinAndMyelinWater_histCompartmentAndCrossR1",'-struct' ...
        ,'r1Data');
end
%% functions in use
function plt = plotVerticalLine(xPos,yValue)
plt = plot([xPos xPos],[0 yValue],'--k','LineWidth',0.7);
end









