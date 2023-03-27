clc; clear all; close all;

%% paths
resultsPath = ...
    "C:\Users\maxoe\Google Drive\Promotion\Simulation\RESULTS\myelinModelResults\";
compositionsFileName = "20230323_whiteMatterCompositionBasedOnLiterature";
compartmentR1FileName = "20230324_compartmentAndCrossR1FromSimulations";

%% load data
compositions = load(resultsPath + compositionsFileName + ".mat");
compartmentR1 = load(resultsPath + compartmentR1FileName + ".mat");

%% litature inputs

% ---- R1 values for the free water part -> CSF
% based on Tsukiashi2018
% T1 vs. temperature: y = 0.058071*x + 1.7623 with x = Temperature in °C
temperature = 36.6; % °C - based on Obermeyer2017
freeWaterR1 =  1/(0.058071 * temperature + 1.7623);
fprintf("Free water R1: %.4f \n", freeWaterR1);
freeWaterR1 = 1/2.5; %CSF value for R1
freeWaterR1 = 0.35;
% 
% % ---- What Heiko calculated
% f_free = 0.8;
% alpha = 0.4;
% r1Lip = 3.

% ---- R1 values measured:
fieldStrengthsInLiterature = ["0.35" "0.55" "1.5" "3" "7"];
r1_WM = cell(length(fieldStrengthsInLiterature),2);
r1_GM = cell(length(fieldStrengthsInLiterature),2);
counter = 1;
for fieldStrength = fieldStrengthsInLiterature
    r1_WM{counter} = fieldStrength;
    r1_GM{counter} = fieldStrength;
    counter = counter + 1;
end

% Deoni2004 - 1.5T:
r1_WM = setR1InDictionary(r1_WM,"1.5",1/0.608);
r1_GM = setR1InDictionary(r1_GM,"1.5",1/1.065);

% Schyboll2018 - 3T
r1_WM = setR1InDictionary(r1_WM,"3",1/0.8);
r1_GM = setR1InDictionary(r1_GM,"3",1/1.25);

% Wang2020 - 0.55T, 1.5T, 3T and 7T
r1_WM = setR1InDictionary(r1_WM,"0.55",1/0.3);
r1_GM = setR1InDictionary(r1_GM,"0.55",1/0.6);

r1_WM = setR1InDictionary(r1_WM,"1.5",1/0.6);
r1_GM = setR1InDictionary(r1_GM,"1.5",1/1);

r1_WM = setR1InDictionary(r1_WM,"3",1/0.85);
r1_GM = setR1InDictionary(r1_GM,"3",1/1.25);

r1_WM = setR1InDictionary(r1_WM,"7",1/1.25);
r1_GM = setR1InDictionary(r1_GM,"7",1/2);


%% Prediction of protein pool R1
fieldStrengthsInData = string(compartmentR1.fieldStrengths);


for counter = 1:length(fieldStrengthsInLiterature)
    fprintf("\n");
    
    fieldStrength = fieldStrengthsInLiterature(counter);
    fprintf("Field strength: %s T.\n",fieldStrength);
    
    r1WmObserved = mean(r1_WM{[r1_WM{:,1}] == fieldStrength,2});
    r1GmObserved = mean(r1_GM{[r1_WM{:,1}] == fieldStrength,2});
    if isnan(r1WmObserved) || isnan(r1GmObserved)
        fprintf("  No literature data for this field strength.\n");
        continue;
    end
    
    fieldStrengthIndex = find(fieldStrengthsInData == fieldStrength);
    if isempty(fieldStrengthIndex)
        fprintf("  No simulated data for this field strength.\n");
        continue;
    end
    
    r1SmSim = compartmentR1.r1_SM(fieldStrengthIndex) ...
        + compartmentR1.r1Auto_SM(fieldStrengthIndex);
    
    % WM
    freeWaterContentWm = compositions.wmWaterContent;
    lipidContentWm = compositions.wmLipidContent;
    r1WmProt = (r1WmObserved - freeWaterContentWm*freeWaterR1) ...
        /((1-freeWaterContentWm)*(1-lipidContentWm)) - lipidContentWm ...
        /(1-lipidContentWm)*r1SmSim;
    fprintf("R1_prot_WM: %.4f\n",r1WmProt);
    
    % GM
    freeWaterContentGm = compositions.gmWaterContent;
    lipidContentGm = compositions.gmLipidContent;
    r1GmProt = (r1GmObserved - freeWaterContentGm*freeWaterR1) ...
        /((1-freeWaterContentGm)*(1-lipidContentGm)) - lipidContentGm ...
        /(1-lipidContentGm)*r1SmSim;
    fprintf("R1_prot_GM: %.4f\n",r1GmProt);
    
end


%% white matter

% T3Index = find(fieldStrengths > 2.999 & fieldStrengths < 3.001);
% R1_obs = 1/0.8;
% f_free = compositions.wmWaterContent;
% alpha = compositions.wmLipidContent;
% r1_lipid = compartmentR1.r1_SM(T3Index) + compartmentR1.r1Auto_SM(T3Index);

% r1_prot = (R1_obs - f_free*freeWaterR1)/((1-f_free)*(1-alpha)) ...
%     - alpha/(1-alpha)*r1_lipid
% 
% 
% %% Gray matter
% sollR1_prot = 2.3;
% 
% fieldStrengths = compartmentR1.fieldStrengths;
% T3Index = find(fieldStrengths > 2.999 & fieldStrengths < 3.001);
% R1_obs = 1/1.2;
% f_free = 0.8;
% alpha = 0.4;
% r1_lipid = compartmentR1.r1_SM(T3Index) + compartmentR1.r1Auto_SM(T3Index);
% 
% r1_prot = (R1_obs - f_free*freeWaterR1)/((1-f_free)*(1-alpha)) ...
%     - alpha/(1-alpha)*r1_lipid


%% functions
function r1Dictionary = setR1InDictionary(r1Dictionary,fieldStrength,value)

r1Dictionary{[r1Dictionary{:,1}] == string(fieldStrength),2}(end+1) ...
    = value;

end