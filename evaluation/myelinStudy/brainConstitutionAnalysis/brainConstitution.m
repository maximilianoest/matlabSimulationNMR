clc; clear all; close all;

%% paths
resultsPath = ...
    "C:\Users\maxoe\Google Drive\Promotion\Simulation\RESULTS" ...
    + "\wholeMyelin_brainConstitution\";

%% data collection

% Data from OBrien1965
% measured at age 10 monthts, 6,9,55 years
data(1).citationKey = "OBrien1965";

data(1).wmLipidContent = [0.490 0.584 0.663 0.646]; 
data(1).wmLipidContent = [0.663 0.646]; % avoid young child data
% data(1).wmLipidContent = [];
data(1).wmProteinContent = 1 - data(1).wmLipidContent;
data(1).wmWaterContent = [0.808 0.755 0.774 0.752]; 
data(1).wmWaterContent = [0.752]; % avoid child data
% data(1).wmWaterContent = [];
data(1).wmNonWaterContent = 1 - data(1).wmWaterContent;

data(1).gmLipidContent = [0.364 0.358 0.376 0.396];
data(1).gmProteinContent = 1 - data(1).gmLipidContent;
data(1).gmWaterContent = [0.841 0.832 0.858 0.823];
data(1).gmNonWaterContent = 1 - data(1).gmWaterContent;

data(1).smLipidContent = [0.780 0.809 0.780 0.780];
data(1).smProteinContent = 1 - data(1).smLipidContent;

% Data from Randall1938
% measured in the following regions:
data(2).citationKey = "Randall1938";
% measurd chemical constituents of twenty three brains from normal and
% psychotic subjects.
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
data(2).wmLipidContent = repmat([0.5730 0.5537 0.5718],1,23);
data(2).wmProteinContent = 1 - data(2).wmLipidContent;
data(2).wmWaterContent = repmat([0.6983 0.7064 0.6992],1,23);
data(2).wmNonWaterContent = 1 - data(2).wmWaterContent;

% [caudate nucleus, frontal cortex, parietal cortex ]
data(2).gmLipidContent = repmat([0.3346 0.3198 0.3277],1,23);
data(2).gmProteinContent = 1 - data(2).gmLipidContent;
data(2).gmWaterContent = repmat([0.8143 0.8412 0.8346],1,23);
data(2).gmNonWaterContent = 1 - data(2).gmWaterContent;


% Data from DeVries1981
% data for WM ( = axolemma-enriched) and myelin exist but not for GM
data(3).citationKey = "DeVries1981";
% density gradient to obtain membrane fraction at the interface of 0.8/1.0
% and 1.0/1.2
% data(3).wmLipidContent = [0.570 0.473];
% data(3).wmProteinContent = 1 - data(3).wmLipidContent;

data(3).smLipidContent = repmat(0.689,1,3);
data(3).smProteinContent = 1 - data(3).smLipidContent;

% Data from Norton:
% The data points from their studies are not included here because they
% investigate bovine brain WM composition.


% Data from Drenthen2021: reference data for training an artifical neural
% network. Data in Table 1 and 2 and REF. 13 volunteers were measured and
% from each volunteer 26 slices were acquiered. REF dataset were 5 slices
% of these 26. Thus, the factor is 13*5 = 65. Given is MWF -> myelin water
% T2 component (15-40ms) divided by all T2 components (15ms - 2000ms) =>
% myelin water divided by whole water == myelin water
data(4).citationKey = "Derenthen2021";
data(4).myelinWaterContent = repmat(0.118,1,13*5);

% Data from Laule2004: measured MWF data from healthy and MS patients. 5 WM
% structures: genu and splenium of CC, posterior internal capsules, minor
% forceps, major forceps. 18 controls -> normal WM MWF
data(5).citationKey = "Laule2004";
data(5).myelinWaterContent = repmat([0.101 0.133 0.156 0.073 0.095],1,18);

%% Data analysis
fieldNames = string(fieldnames(data))';
data(end+1).citationKey = "Averaged";
allData = nan(length(fieldNames)-1,200);
fieldNameCounter = 1;

initializeFigure();
lgd = legend();
lgd.Visible = 'off';

fprintf("<strong>On average </strong> \n");
for fieldName = fieldNames
    if fieldName == "citationKey"
        continue;
    end
    dataArray = [data(:).(fieldName)];
   
    data(end).(fieldName) = mean(dataArray);
    averageData.(fieldName) = data(end).(fieldName);
    fprintf("%-20s: %.4f \n\n",fieldName,data(end).(fieldName));
    
    allData(fieldNameCounter,1:length(dataArray)) = dataArray;
    fieldNameCounter = fieldNameCounter + 1;
end

boxplot(allData','Labels',fieldNames(2:end));
xtickangle(45)


save(resultsPath + datestr(now,"yyyymmdd") + "_" ...
    + "wmAndGMCompositionBasedOnLiterature",'-struct','averageData');




















