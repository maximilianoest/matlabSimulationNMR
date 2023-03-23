clc; clear all; close all;

%% paths
resultsPath = ...
    "C:\Users\maxoe\Google Drive\Promotion\Simulation\RESULTS\myelinModelResults\";

%% data collection

% Data from OBrien1965
% measured at age 10 monthts, 6,9,55 years
data(1).citationKey = "OBrien1965";

data(1).wmLipidContent = [0.490 0.584 0.663 0.646]; 
data(1).wmLipidContent = [0.663 0.646]; % avoid child data
data(1).wmNonLipidContent = 1 - data(1).wmLipidContent;
data(1).wmWaterContent = [0.808 0.755 0.774 0.752]; 
data(1).wmWaterContent = [0.774 0.752]; % avoid child data
data(1).wmNonWaterContent = 1 - data(1).wmWaterContent;

data(1).gmLipidContent = [0.364 0.358 0.376 0.396];
data(1).gmNonLipidContent = 1 - data(1).gmLipidContent;
data(1).gmWaterContent = [0.841 0.832 0.858 0.823];
data(1).gmNonWaterContent = 1 - data(1).gmWaterContent;

data(1).smLipidContent = [0.780 0.809 0.780 0.780];
data(1).smNonLipidContent = 1 - data(1).smLipidContent;

% Data from Randall1938
% measured in the following regions:
data(2).citationKey = "Randall1938";
% WM
% - corona radiata = white matter containing the great conduction path ways
% - frontal white matter = pathways from frontal cortex cell bodies to
% cortex and basal ganglia
% - parietal white matter = connecting subcortical regions with parietal,
% occipital and temporal cortex
% GM & WM
% - brain stem: below main brain without cerebellum -> contains WM and GM 
% structures -> resembles white matter with significant differences
% - thalamus: before nerve fibers from senses enter the cortex they run
% through the thalamus -> contains WM and GM structures -> resembles gray
% matter with significant difference
% GM
% - caudate nucleus: regions with nuclei as part of the basal ganglia ->
% gray matter
% - frontal cortex: mainly gray nerve cells on top of frontal WM -> GM
% - parietal cortex: mainly gray nerve cells on top of parietal WM -> GM

% [corona radiata, frontal white, parietal white]
data(2).wmLipidContent = [0.5730 0.5537 0.5718];
data(2).wmNonLipidContent = 1 - data(2).wmLipidContent;
data(2).wmWaterContent = [0.6983 0.7064 0.6992];
data(2).wmNonWaterContent = 1 - data(2).wmWaterContent;

% [caudate nucleus, frontal cortex, parietal cortex
data(2).gmLipidContent = [0.3346 0.3198 0.3277];
data(2).gmNonLipidContent = 1 - data(2).gmLipidContent;
data(2).gmWaterContent = [0.8143 0.8412 0.8346];
data(2).gmNonWaterContent = 1 - data(2).gmWaterContent;


% Data from DeVries1981
% data for WM ( = axolemma-enriched) and myelin exist but not for GM
data(3).citationKey = "DeVries1981";
% density gradient to obtain membrane fraction at the interface of 0.8/1.0
% and 1.0/1.2
data(3).wmLipidContent = [0.570 0.473];
data(3).wmNonLipidContent = 1 - data(3).wmLipidContent;

data(3).smLipidContent = 0.689;
data(3).smNonLipidContent = 1 - data(3).smLipidContent;

% Data from Norton:
% The data points from their studies are not included here because they
% investigate bovine brain WM composition.


%% Data analysis
fieldNames = string(fieldnames(data))';
data(end+1).citationKey = "Averaged";

fprintf("<strong>On average </strong> \n");
for fieldName = fieldNames
    if fieldName == "citationKey"
        continue;
    end
   data(end).(fieldName) = mean([data(:).(fieldName)]);
   averageData.(fieldName) = data(end).(fieldName);
   fprintf("%-20s: %.4f \n\n",fieldName, data(end).(fieldName))
end


save(resultsPath + "_" + datestr(now,"yyyymmdd") ...
    + "whiteMatterCompositionBasedOnLiterature",'-struct','averageData');




















