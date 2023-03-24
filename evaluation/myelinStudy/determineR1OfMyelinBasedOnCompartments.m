clc; clear all; close all;

clc; clear all; close all; fclose('all');
addpath(genpath("../../library/"));
constants = readConstantsFile(sprintf("..%s..%stxtFiles%sconstants.txt" ...
    ,createFilesepStringArray(3)));


resultsDir = "C:\Users\maxoe\Google Drive\Promotion\Simulation\RESULTS\";
folderNames = ["solidMyelin_20230215_determineShortPaddedCorrFunc" ...
    "myelinWater_20230206_DetermineCorrelationFunctions" ...
    "wholeMyelin_20230302_determineCrossRelaxationCorrFuncs"] + filesep;
fileNames = ["20230215_Results_MYELINsolidMyelin_20221222_MYELIN_TIP4_Bilayer_50water_solidMyelin_H_whole_dt4ps_simTime800ns" ...
    "20230208_Results_MYELINmyelinWater_20230202_MYELIN_TIP4_Monolayer_50water_myelinWater_H_whole_dt01ps_simTime10ns" ...
    "20230310_Results_MYELINallMonolayer_20230202_MYELIN_TIP4_Monolayer_50water_allMonolayer_H_whole_dt01ps_simTime5ns"] ...
    + ".mat";
compartmentNames = ["solid" "water" "crossAuto"];

filePaths = checkIfFilesExistAndCreateFilePaths(resultsDir,folderNames ...
    ,fileNames);
dataArray = loadFiles(filePaths);

savingDir = resultsDir + "myelinModelResults\";
saving = 1;

%% ---- configuration
dipolDipolConstant = 3/4*(constants.vaccumPermeability/(4*pi) ...
    *constants.hbar*constants.gyromagneticRatioOfHydrogenAtom^2)^2 ...
    /(constants.nanoMeter^6);
fieldStrengths = 0.1 : 0.05 : 7;
larmorFrequencies = constants.gyromagneticRatioOfHydrogenAtom ...
    *fieldStrengths;


%% ---- solid myelin
data = dataArray{2,[dataArray{1,:}] == "solidMyelin"};
fprintf("<strong> --- %s : </strong>\n",data.constituent);
offsetCorrectionRegion = [450e-9 550e-9];
cutOffTime_SM = 600e-9;

nNCases_SM = data.nearestNeighbourCases;
[~,nNIdx] = max(nNCases_SM);
fprintf("Nearest neighbor cases: [ %s] \n",sprintf("%i ",nNCases_SM));
deltaTInS_SM = data.deltaTInS;
fprintf("dT: %.4d sec \n",deltaTInS_SM);
indices = offsetCorrectionRegion/deltaTInS_SM;
cutOffIndex_SM = cutOffTime_SM/deltaTInS_SM;

corrFuncFirstOrder_SM = correctOffsetInCorrFunc( ...
    data.sumCorrFuncFirstOrder(nNIdx,:)/data.atomCounter,indices);
corrFuncFirstOrder_SM(cutOffIndex_SM+1:end) = 0;
corrFuncSecondOrder_SM = correctOffsetInCorrFunc( ...
    data.sumCorrFuncSecondOrder(nNIdx,:)/data.atomCounter,indices);
corrFuncSecondOrder_SM(cutOffIndex_SM+1:end) = 0;


%% ---- myelin water
data = dataArray{2,[dataArray{1,:}] == "myelinWater"};
fprintf("<strong> --- %s : </strong>\n",data.constituent);
cutOffTime_MW = 6e-9;

nNCases_MW = data.nearestNeighbourCases;
[~,nNIdx] = max(nNCases_MW);
fprintf("Nearest neighbor cases: [ %s] \n",sprintf("%i ",nNCases_MW));
deltaTInS_MW = data.deltaTInS;
fprintf("dT: %.4d sec \n",deltaTInS_MW);
cutOffIndex_MW = cutOffTime_MW/deltaTInS_MW;

corrFuncFirstOrder_MW = data.sumCorrFuncFirstOrder(nNIdx,:) ...
    /data.atomCounter;
corrFuncFirstOrder_MW(cutOffIndex_MW+1:end) = 0;
corrFuncSecondOrder_MW = data.sumCorrFuncSecondOrder(nNIdx,:) ...
    /data.atomCounter;
corrFuncSecondOrder_MW(cutOffIndex_MW+1:end) = 0;

%% ---- cross and auto relaxation
data = dataArray{2,[dataArray{1,:}] == "allMonolayer"};
fprintf("<strong> --- %s : </strong>\n",data.constituent);

nNCases_CA = data.nearestNeighbourCases;
[~,nNIdx] = max(nNCases_CA);
fprintf("Nearest neighbor cases: [ %s] \n",sprintf("%i ",nNCases_CA));
deltaTInS_CA = data.deltaTInS;
fprintf("dT: %.4d sec \n",deltaTInS_CA);

corrFuncZerothOrder_CA = data.sumCorrFuncZerothOrder(nNIdx,:) ...
    /data.atomCounter;
corrFuncFirstOrder_CA = data.sumCorrFuncFirstOrder(nNIdx,:) ...
    /data.atomCounter;
corrFuncSecondOrder_CA = data.sumCorrFuncSecondOrder(nNIdx,:) ...
    /data.atomCounter;

hCountInMembrane = length(data.matFileIndices{ ...
    data.groupsToSearchIn == "MEMB"});
hCountInWater = length(data.matFileIndices{ ...
    data.groupsToSearchIn == "Water"});


%% ---- calaculate R1
r1_SM = zeros(1,length(fieldStrengths));
r1_MW = zeros(1,length(fieldStrengths));

r1Auto_SM = zeros(1,length(fieldStrengths));
r1Cross_SM = zeros(1,length(fieldStrengths));

for fieldStrengthNr = 1:length(fieldStrengths)
    % solid myelin
    [specDensFirstOrder,specDensSecondOrder] = ...
        calculateSpectralDensities(corrFuncFirstOrder_SM ...
        ,corrFuncSecondOrder_SM,larmorFrequencies(fieldStrengthNr) ...
        ,deltaTInS_SM,[0.9 1]);
    r1_SM(fieldStrengthNr) = calculateR1WithSpectralDensity( ...
        specDensFirstOrder,specDensSecondOrder,dipolDipolConstant);
    
    % myelin water
    [specDensFirstOrder,specDensSecondOrder] = ...
        calculateSpectralDensities(corrFuncFirstOrder_MW ...
        ,corrFuncSecondOrder_MW,larmorFrequencies(fieldStrengthNr) ...
        ,deltaTInS_MW,[0.9 1]);
    r1_MW(fieldStrengthNr) = calculateR1WithSpectralDensity( ...
        specDensFirstOrder,specDensSecondOrder,dipolDipolConstant);
    
    % auto relaxtion
    [specDensZerothsOrder,specDensFirstOrder,specDensSecondOrder] = ...
        calculateSpecDensForZerothFirstAndSecondOrder( ...
        corrFuncZerothOrder_CA,corrFuncFirstOrder_CA ...
        ,corrFuncSecondOrder_CA,larmorFrequencies(fieldStrengthNr) ...
        ,deltaTInS_CA,[0.9 1]);
    
    r1Auto_SM(fieldStrengthNr) = calculateR1AutoWithSpecDens( ...
        specDensZerothsOrder,specDensFirstOrder,specDensSecondOrder ...
        ,dipolDipolConstant);
    
    % cross relaxation 
    r1Cross_SM(fieldStrengthNr) = calculateR1CrossWithSpecDens( ...
        specDensZerothsOrder,specDensSecondOrder,dipolDipolConstant);
    
    
    
    
    
end

% effective relaxation rates
r1Eff_SM = r1_SM + r1Auto_SM;

r1Auto_MW = hCountInMembrane/hCountInWater * r1Auto_SM;
r1Eff_MW = r1_MW + r1Auto_MW;
r1Cross_MW = hCountInMembrane/hCountInWater * r1Cross_SM;
%% ---- plotting
fig1 = initializeFigure();
subPlt1 = initializeSubplot(fig1,2,2,1);
subPlt2 = initializeSubplot(fig1,2,2,2);
subPlt3 = initializeSubplot(fig1,2,2,3);
subPlt4 = initializeSubplot(fig1,2,2,4);

%% ---- solid myelin
fig1.CurrentAxes = subPlt1;
timeAxis_SM = 0:deltaTInS_SM:(length(corrFuncFirstOrder_SM)-1)*deltaTInS_SM;
plot(timeAxis_SM,real(corrFuncFirstOrder_SM));
plot(timeAxis_SM,imag(corrFuncFirstOrder_SM));
plot(timeAxis_SM,real(corrFuncSecondOrder_SM));
plot(timeAxis_SM,imag(corrFuncSecondOrder_SM));

axis([0 cutOffTime_SM -inf 4e3]);
legend("real 1$^{st}$","imag 1$^{st}$","real 2$^{nd}$","imag 2$^{nd}$");
title("Corr. func. SM");
xlabel("Delay time [s]");

%% ---- myelin water
fig1.CurrentAxes = subPlt2;
timeAxis_MW = 0:deltaTInS_MW:(length(corrFuncFirstOrder_MW)-1)*deltaTInS_MW;
plot(timeAxis_MW,real(corrFuncFirstOrder_MW));
plot(timeAxis_MW,imag(corrFuncFirstOrder_MW));
plot(timeAxis_MW,real(corrFuncSecondOrder_MW));
plot(timeAxis_MW,imag(corrFuncSecondOrder_MW));

axis([0 cutOffTime_MW -inf 1.5e3]);
legend("real 1$^{st}$","imag 1$^{st}$","real 2$^{nd}$","imag 2$^{nd}$");
title("Corr. func. MW");
xlabel("Delay time [s]")


%% ---- correlation functions for R1_SM^cross and ^auto
fig1.CurrentAxes = subPlt3;
timeAxis_CA = 0:deltaTInS_CA:(length(corrFuncFirstOrder_CA)-1)*deltaTInS_CA;
plot(timeAxis_CA,real(corrFuncZerothOrder_CA));
plot(timeAxis_CA,imag(corrFuncZerothOrder_CA));
plot(timeAxis_CA,real(corrFuncFirstOrder_CA));
plot(timeAxis_CA,imag(corrFuncFirstOrder_CA));
plot(timeAxis_CA,real(corrFuncSecondOrder_CA));
plot(timeAxis_CA,imag(corrFuncSecondOrder_CA));

legend("real 0$^{th}$","imag 0$^{th}$","real 1$^{st}$" ...
    ,"imag 1$^{st}$","real 2$^{nd}$","imag 2$^{nd}$");
title("Corr. func. for cross and auto");
xlabel("Delay time [s]")

%% ---- r1 rates
fig1.CurrentAxes = subPlt4;
plot(fieldStrengths,r1Eff_SM);
plot(fieldStrengths,r1Eff_MW);

plot(fieldStrengths,r1_SM);
plot(fieldStrengths,r1_MW);

plot(fieldStrengths,r1Auto_SM,'--');
plot(fieldStrengths,r1Cross_SM,':');

plot(fieldStrengths,r1Auto_MW,'--');
plot(fieldStrengths,r1Cross_MW,':');

axis([0 inf 0 10]);
legend("R$^{eff}_{1,SM}$","R$^{eff}_{1,MW}$","R$_{1,SM}$","R$_{1,MW}$" ...
    ,"R$^{auto}_{1,SM}$","R$^{cross}_{1,SM}$","R$^{auto}_{1,MW}$" ...
    ,"R$^{cross}_{1,MW}$");
xticks(0:0.5:7);
grid on

fig2 = initializeFigure();
plot(fieldStrengths,r1Eff_SM);
plot(fieldStrengths,r1Eff_MW);

plot(fieldStrengths,r1_SM);
plot(fieldStrengths,r1_MW);

plot(fieldStrengths,r1Auto_SM,'--');
plot(fieldStrengths,r1Cross_SM,':');

plot(fieldStrengths,r1Auto_MW,'--');
plot(fieldStrengths,r1Cross_MW,':');

axis([0 inf 0 10]);
legend("R$^{eff}_{1,SM}$","R$^{eff}_{1,MW}$","R$_{1,SM}$","R$_{1,MW}$" ...
    ,"R$^{auto}_{1,SM}$","R$^{cross}_{1,SM}$","R$^{auto}_{1,MW}$" ...
    ,"R$^{cross}_{1,MW}$");
xticks(0:0.5:7);


%% ---- saving

if saving
    set(0,'CurrentFigure',fig1);
    saveFigureTo(savingDir,"wholeMyelin",datestr(now,'yyyymmdd') ...
        ,"CompartmentAndCrossR1",true);
end

save(resultsDir + "myelinModelResults\" + datestr(now,"yyyymmdd") ...
    +"_compartmentAndCrossR1FromSimulations",'r1_MW','r1_SM' ...
    ,'r1Auto_MW','r1Auto_SM','r1Cross_MW','r1Cross_SM','r1Eff_MW' ...
    ,'r1Eff_SM','fieldStrengths');





%% ---- functions
function filePaths = checkIfFilesExistAndCreateFilePaths(resultsDir ...
    ,folderNames,fileNames)
filePaths = strings(1,length(fileNames));
for fileNr = 1:length(fileNames)
    filePath = resultsDir + folderNames(fileNr) + fileNames(fileNr);
    if exist(filePath,'file')
        filePaths(fileNr) = filePath;
    else
        error("One file does not exist! %s",filePath);
    end
end

end

function data = loadFiles(filePaths)

data = cell(2,length(filePaths));

for fileNr = 1:length(filePaths)
    
    data{2,fileNr} = load(filePaths(fileNr));
    data{1,fileNr} = string(data{2,fileNr}.constituent);
end

end

function offsetCorrectedCorrFunc = correctOffsetInCorrFunc(corrFunc ...
    ,indices)
offsetCorrectedCorrFunc = corrFunc - squeeze(mean(corrFunc( ...
    indices(1):indices(2))));

end


