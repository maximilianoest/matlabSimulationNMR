clc; clear all; close all; fclose('all');

%% set dependencies
addpath(genpath(sprintf('..%s..%s..%slibrary',filesep,filesep,filesep)));
addpath(genpath(sprintf('..%s..%s..%stxtFiles' ...
    ,createFilesepStringArray(3))));
constants = readConstantsFile('constants.txt');

%% set up
saving = 0;
smallFieldStrength = 0.6; % T

%% directories
resultsDir = sprintf("..%s..%s..%sRESULTS%",createFilesepStringArray(3));

%% constitutions
constitutionFolderName = sprintf("%swholeMyelin_brainConstitution%s" ...
    ,createFilesepStringArray(2));
constitutionFileName = "20230419_wmAndGMCompositionBasedOnLiterature";
constData = load(resultsDir + constitutionFolderName ...
    + constitutionFileName);

%% relaxation rates 
r1RatesFolder = sprintf("%swholeMyelin_relaxationRates%s" ...
    ,createFilesepStringArray(2));
r1RatesFileName  ...
    = "20230419_solidMyelinAndMyelinWater_histCompartmentAndCrossR1";
r1Data =  load(resultsDir + r1RatesFolder + r1RatesFileName);

%% observable relaxation rates in qMRI
tissueR1RatesFileName = "20230419_tissueR1BasedOnLiterature";
tissueR1Data = load(resultsDir + r1RatesFolder + tissueR1RatesFileName);

[wmFieldStrengthsInLiterature,relaxationRatesWM] = ...
    getFieldStrengthsAndCalculateAvgR1(tissueR1Data.r1_WM);
[gmFieldStrengthsInLiterature,relaxationRatesGM] = ...
    getFieldStrengthsAndCalculateAvgR1(tissueR1Data.r1_GM);

freeWaterR1 = r1Data.r1_free;

%% calculate protein pool R1 based on observations and constitutions

% Based on the assumptions:
% -fast exchange between all components => observed R1 can be seen as
% linear combination of compartment R1. The coefficients are given by the
% fraction of each compartment
% - solid myelin R1 is equal for GM and WM 
% - identical surface water R1 in protein and lipid


if ~any(gmFieldStrengthsInLiterature == wmFieldStrengthsInLiterature)
    error("Not implemented case for literature data");
end

fieldStrengthsInLit = wmFieldStrengthsInLiterature;
simFieldStrengths = r1Data.fieldStrengths;
r1Data.litFieldStrengths = zeros(1,length(fieldStrengthsInLit));
r1Data.r1Observed_WM = zeros(1,length(fieldStrengthsInLit));
r1Data.r1Observed_GM = zeros(1,length(fieldStrengthsInLit));

r1Data.r1Prot_WM = zeros(1,length(fieldStrengthsInLit));
r1Data.r1Prot_GM = zeros(1,length(fieldStrengthsInLit));

% WM fractions
wmLipidWaterContent = constData.wmWaterContent ...
    * constData.myelinWaterContent;
wmSolidLipidContent = constData.wmNonWaterContent ...
    * constData.wmLipidContent;
wmESIVR = wmLipidWaterContent/wmSolidLipidContent;
wmSolidProteinContent = constData.wmNonWaterContent ...
    * constData.wmProteinContent;
wmProtWaterContent = wmSolidProteinContent*wmESIVR;
wmFreeWaterContent = constData.wmWaterContent - wmLipidWaterContent ...
    - wmProtWaterContent;

% GM fractions
gmLipidWaterContent = constData.gmNonWaterContent ...
    * constData.gmLipidContent * wmESIVR;
gmSolidLipidContent = constData.gmNonWaterContent ...
    * constData.gmLipidContent;
gmProtWaterContent = constData.gmNonWaterContent ...
    * constData.gmProteinContent * wmESIVR;
gmSolidProteinContent = constData.gmNonWaterContent ...
    * constData.gmProteinContent;
gmFreeWaterContent = constData.gmWaterContent - gmLipidWaterContent ...
    - gmProtWaterContent;



% WM for small field strengths where lipid cut out
wmNonLipidContent = 1 - wmSolidLipidContent;
smallFieldWmSolidProteinContent = wmSolidProteinContent  ...
    / wmNonLipidContent;
smallFieldWmLipidWaterContent = wmLipidWaterContent / wmNonLipidContent;
smallFieldWmProteinWaterContent = wmProtWaterContent / wmNonLipidContent;
smallFieldFreeWaterContent = wmFreeWaterContent / wmNonLipidContent;


% GM for small field strengths where lipid is cut out
gmNonLipidContent = 1 - wmSolidLipidContent;
smallFieldGmSolidProteinContent = gmSolidProteinContent ...
    / gmNonLipidContent;
smallFieldGmLipidWaterContent = gmLipidWaterContent / gmNonLipidContent;
smallFieldGmProteinWaterContent = gmProtWaterContent / gmNonLipidContent;
smallFieldGmFreeWaterContent = gmFreeWaterContent / gmNonLipidContent;


for fieldStrengthNr = 1:length(fieldStrengthsInLit)
    fieldStrength = fieldStrengthsInLit(fieldStrengthNr);
    fprintf("<strong>Field strength: %.4f T\n</strong>",fieldStrength);
    r1Data.litFieldStrengths(fieldStrengthNr) = fieldStrength;
    
    r1_LW = r1Data.r1Hist_MW(simFieldStrengths > (fieldStrength-0.001) ...
        & simFieldStrengths < (fieldStrength+0.001));
    
    r1_SL = r1Data.r1Eff_SM(simFieldStrengths > (fieldStrength-0.001) ...
        & simFieldStrengths < (fieldStrength+0.001));
    fprintf("  SL R1 MD-Sim: %.4f\n",r1_SL);
    
    if length(r1_LW) ~= 1 || length(r1_SL) ~= 1
        error("More than 2 R1 are found for a field strength.");
    end
    
    % WM protein R1
    observedWMR1 = relaxationRatesWM(fieldStrengthNr);
    r1Data.r1Observed_WM(fieldStrengthNr) = observedWMR1;
    
    % GM protein R1
    observedGMR1 = relaxationRatesGM(fieldStrengthNr);
    r1Data.r1Observed_GM(fieldStrengthNr) = observedGMR1;
    
    
    if fieldStrength < smallFieldStrength
        % In the low field case, the lipid pool is cut out. The field
        % strength under which the lipid pool is cut out is determined
        % under set up.
        
        
        r1Data.smallFieldStrengths = fieldStrength;
        r1Data.r1Prot_WM(fieldStrengthNr) = ...
            1/(smallFieldWmSolidProteinContent) ...
            *(observedWMR1 ...
            - smallFieldWmLipidWaterContent*r1_LW ...
            - smallFieldWmProteinWaterContent*r1_LW ...
            - smallFieldFreeWaterContent*freeWaterR1);
        fprintf("  Protein R1 WM (without lipid): %.4f. \n" ...
            ,r1Data.r1Prot_WM(fieldStrengthNr));
        
        r1Data.r1Prot_GM(fieldStrengthNr) = ...
            1/(smallFieldGmSolidProteinContent) ...
            *(observedGMR1 ...
            - smallFieldGmLipidWaterContent*r1_LW ...
            - smallFieldGmProteinWaterContent*r1_LW ...
            - smallFieldGmFreeWaterContent*freeWaterR1);
        fprintf("  Protein R1 GM (without lipid): %.4f. \n" ...
            ,r1Data.r1Prot_GM(fieldStrengthNr));
    else
        r1Data.r1Prot_WM(fieldStrengthNr) = ...
            1/wmSolidProteinContent ...
            *(observedWMR1 ...
            - wmSolidLipidContent*r1_SL...
            - wmLipidWaterContent*r1_LW ...
            - wmProtWaterContent*r1_LW ...
            - wmFreeWaterContent*freeWaterR1);
        fprintf("  Protein R1 WM: %.4f. \n",r1Data.r1Prot_WM( ...
            fieldStrengthNr));
        
        r1Data.r1Prot_GM(fieldStrengthNr) = ...
            1/gmSolidProteinContent ...
            *(observedGMR1 ...
            - gmSolidLipidContent*r1_SL ...
            - gmLipidWaterContent*r1_LW ...
            - gmProtWaterContent*r1_LW ...
            - gmFreeWaterContent*freeWaterR1);
        fprintf("  Protein R1 GM: %.4f. \n",r1Data.r1Prot_GM( ...
            fieldStrengthNr));
    end
    
    
end

if saving
    save(resultsDir +r1RatesFolder + datestr(now,"yyyymmdd") ...
        + "_SM_MW_SP_histCompartmentAndCrossR1",'-struct','r1Data');
end


%% functions
function [fieldStrengths,r1Rates] = getFieldStrengthsAndCalculateAvgR1( ...
    cellArray)
fieldStrengths = zeros(1,size(cellArray,1));
r1Rates = zeros(1,size(cellArray,1));

for lineNr = 1:size(cellArray,1)
    
    fieldStrengths(lineNr) = str2double(cellArray{lineNr,1});
    if isnan(fieldStrengths(lineNr))
       warning("The field strength %s cannot be converted" ...
           ,cellArray{lineNr,1});
    end
    
    r1Rates(lineNr) = mean(cellArray{lineNr,2});
    if isempty(r1Rates(lineNr))
        warning("The field strength %s has no relaxation rate value." ...
            ,r1Rates(lineNr));
    end
end

end











