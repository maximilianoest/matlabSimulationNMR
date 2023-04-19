clc; clear all; close all;

%% paths
resultsPath = ...
    "C:\Users\maxoe\Google Drive\Promotion\Simulation\RESULTS\myelinModelResults\";

addpath(genpath(sprintf('..%s..%s..%slibrary',filesep,filesep,filesep)));
addpath(genpath(sprintf('..%s..%s..%stxtFiles' ...
    ,createFilesepStringArray(3))));
constants = readConstantsFile('constants.txt');

resultsDir = sprintf("..%s..%s..%sRESULTS%s" ...
    + "wholeMyelin_relaxationRates%s",createFilesepStringArray(5));

%% litature inputs

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
r1_WM = setR1InDictionary(r1_WM,"0.55",1/0.4);
r1_GM = setR1InDictionary(r1_GM,"0.55",1/0.6);

r1_WM = setR1InDictionary(r1_WM,"1.5",1/0.5);
r1_GM = setR1InDictionary(r1_GM,"1.5",1/0.9);

r1_WM = setR1InDictionary(r1_WM,"3",1/0.85);
r1_GM = setR1InDictionary(r1_GM,"3",1/1.25);

r1_WM = setR1InDictionary(r1_WM,"7",1/1.25);
r1_GM = setR1InDictionary(r1_GM,"7",1/2);

% Heiko's assumptions - 0.35T, 1.5T and 3T
r1_WM = setR1InDictionary(r1_WM,"0.35",1/0.33);
r1_GM = setR1InDictionary(r1_GM,"0.35",1/0.5);

r1_WM = setR1InDictionary(r1_WM,"1.5",1/0.55);
r1_GM = setR1InDictionary(r1_GM,"1.5",1/0.8);

r1_WM = setR1InDictionary(r1_WM,"3",1/0.8);
r1_GM = setR1InDictionary(r1_GM,"3",1/1.2);

save(resultsDir + datestr(now,"yyyymmdd") ...
    + "_tissueR1BasedOnLiterature",'r1_GM','r1_WM');

%% functions
function r1Dictionary = setR1InDictionary(r1Dictionary,fieldStrength,value)

r1Dictionary{[r1Dictionary{:,1}] == string(fieldStrength),2}(end+1) ...
    = value;

end