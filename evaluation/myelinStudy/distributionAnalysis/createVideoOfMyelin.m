clc; clear all; close all; fclose('all');

addpath(genpath(sprintf('..%s..%s..%slibrary',filesep,filesep,filesep)));
addpath(genpath(sprintf('..%s..%s..%stxtFiles' ...
    ,createFilesepStringArray(3))));



resultsDir = sprintf("..%s..%s..%sRESULTS%swholeMyelin" ...
    + "_mdVideo%s",createFilesepStringArray(5));
if ~exist(resultsDir,'dir')
    warning('The results directory does not exist but will be created.');
    mkdir(resultsDir)
end

directory = "C:\Users\maxoe\Documents\Gromacs\testDatasets\";
folderName = "MYELIN";
layerForm = "/Monolayer/";
fileName = "20230404_MYELIN_TIP4_50Water_ShortSimDur_prd_allAtoms_dt0_05ps_simDur0_07046ns";
groFileName = "step6.5_equilibration";

locationData = load(directory + folderName + layerForm + fileName);
meanLocationData = mean(mean(locationData.trajectories,3),1);
locationData.trajectories = locationData.trajectories - meanLocationData;

splittedFileName = strsplit(fileName,'_');
if length(splittedFileName) == 12
    simDurInNs = str2num(strrep(splittedFileName(end-1),'simTime','')) ...
        + str2num(strrep(splittedFileName(end),'ns',''))/1000;
    deltaTInPs = str2num(strrep(splittedFileName(end-3),'dt','')) ...
        + str2num(strrep(splittedFileName(end-2),'ps',''))/1000;
elseif length(splittedFileName) == 11
    simDurInNs = str2num(strrep(splittedFileName(end-1),'simDur','')) ...
        + str2num(strrep(splittedFileName(end),'ns',''))/10000;
    deltaTInPs = str2num(strrep(splittedFileName(end-3),'dt','')) ...
        + str2num(strrep(splittedFileName(end-2),'ps',''))/100;
else
    deltaTInPs = str2num(strrep(strrep(splittedFileName(9),'ps',''),'dt','')); %#ok<ST2NM>
    simDurInNs = str2num(strrep(splittedFileName(10),'simTime','')) ...
        + str2num(strrep(splittedFileName(11),'ns',''))/1000; %#ok<ST2NM>
end

atomNamesToIgnore = ["POT" "CLA" "M"];
moleculesToIgnore = ["POT" "CLA"];

%% --------- find out which index belongs to which atom and molecule
fprintf("Determining which index belongs to which kind of atom. \n");
groFileID = fopen(sprintf("%s%s.gro",directory + folderName ...
    + layerForm,groFileName));
groFileLine = fgetl(groFileID);
splittedLineString = strsplit(groFileLine," ");
if strcmp(splittedLineString{1}, "Title")
    groFileLine = fgetl(groFileID);
    atomsCountFromGroFile = str2num(groFileLine); %#ok<ST2NM>
    atomNameArray = strings(1,atomsCountFromGroFile);
    moleculeNameArray = strings(1,atomsCountFromGroFile);
    groFileLine = fgetl(groFileID);
end

atomNameCollection = strings();
moleculeNameCollection = strings();
oldMoleculeIndex = 1;
oldMoleculeName = string.empty();

moleculeAtomCollection = [];
moleculeStartIndex = 1;
moleculeCounter = zeros(1,11);
moleculeLengths = [];

while groFileLine ~= -1
    splittedLineString = strsplit(groFileLine," ");
    groFileLine = fgetl(groFileID);
    
    splittedLineString = splittedLineString(2:end);
    
    if length(splittedLineString) < 8
        moleculeCounter(moleculeNameCollection == moleculeNameArray( ...
                atomIndex - 1)) = moleculeCounter( ...
                moleculeNameCollection == moleculeNameArray( ...
                atomIndex - 1)) + 1;
        break;
    end
    
    moleculeName = splittedLineString{1};
    
    % check for different line styles
    if length(splittedLineString) == 9
        atomName = splittedLineString{2};
        atomIndex = str2num(splittedLineString{3}); %#ok<ST2NM>
    elseif length(splittedLineString) == 8 && atomIndex >= 9999
        atomAndIndex = splittedLineString{2};
        atomName = atomAndIndex(1:end-5);
        atomIndex = str2double(atomAndIndex(end-4:end));
        if atomIndex > atomsCountFromGroFile
            error(['There is a higher index than there are atoms' ...
                ' in the sample according to gro file.']);
        end
    else
        continue;
    end 
    
    
    % check for different atom descriptions
    if contains(atomName,'POT') || contains(atomName,'CLA')
        if any(atomNamesToIgnore == atomName(1:3))
            atomNameArray(atomIndex) = "ignored";
        else
            atomNameArray(atomIndex) = atomName(1:3);
        end
    else
        if any(atomNamesToIgnore == atomName(1))
            atomNameArray(atomIndex) = "ignored";
        else
            atomNameArray(atomIndex) = atomName(1);
        end
    end
    
    moleculeIndex = getMoleculeIndexFromMoleculeName(moleculeName);
    
    % molecule name collection    
    if contains(moleculeName,"CER24") 
        moleculeName = "GAL";
    elseif contains(moleculeName,"BGAL")
        moleculeName = "BGAL";
    elseif contains(moleculeName,"CHL")
        moleculeName = "CHL";
    elseif contains(moleculeName,"DSPE")
        moleculeName = "DSPE";
    elseif contains(moleculeName,"DOPC")
        moleculeName = "DOPC";
    elseif contains(moleculeName,"SOPS")
        moleculeName = "SOPS";
    elseif contains(moleculeName,"SSM")
        moleculeName = "SSM";
    elseif contains(moleculeName,"PSPI")
        moleculeName = "PSPI";
    elseif contains(moleculeName,"POT")
        moleculeName = "POT";
    elseif contains(moleculeName,"CLA")
        moleculeName = "CLA";
    elseif contains(moleculeName,"SOL")
        moleculeName = "WATER";
    else
        error("Molecule not found");
    end
    
    
    moleculeNameArray(atomIndex) = moleculeName;
    
    if isempty(oldMoleculeName)
        oldMoleculeName = moleculeName;
    end
    
    
    if oldMoleculeIndex ~= moleculeIndex
        moleculeLength = atomIndex - moleculeStartIndex;
        
        if any([1 3] == moleculeIndex) % Then, the molecule = GAL + BGAL
            if moleculeLength == 151
                moleculeNameArray(moleculeStartIndex:atomIndex-1) ...
                    = "GALS";
            else
                moleculeNameArray(moleculeStartIndex:atomIndex-1) ...
                    = "GALC";
            end
        end
        
        % find GALS and GALC
        if moleculeName ~= "BGAL"
            if ~contains(moleculeNameCollection,moleculeNameArray( ...
                    atomIndex - 1))
                if moleculeNameCollection(1) == ""
                    moleculeNameCollection(1) = moleculeNameArray( ...
                        atomIndex - 1);
                else
                    moleculeNameCollection(end+1) = moleculeNameArray( ...
                        atomIndex - 1); %#ok<SAGROW>
                end
            end
            moleculeLengths(end+1) = moleculeLength;
            moleculeCounter(moleculeNameCollection == moleculeNameArray( ...
                atomIndex - 1)) = moleculeCounter( ...
                moleculeNameCollection == moleculeNameArray( ...
                atomIndex - 1)) + 1;
            
            moleculeStartIndex = atomIndex;
            
        end
        oldMoleculeIndex = moleculeIndex;
    end
    
    
    if ~contains(atomNameCollection,atomNameArray(atomIndex))
        if atomNameCollection(1) == ""
            atomNameCollection(end) = atomNameArray(atomIndex);
        else
            atomNameCollection(end+1) = atomNameArray(atomIndex); %#ok<SAGROW>
        end
    end
    
end

%% evaluation

hIndices = atomNameArray == "H";
hLocations = locationData.trajectories(hIndices,:,:);
hColor = [77 85 240]/255; % blue
hSymbole = '*';
hLineWidth = 3;

oIndices = atomNameArray == "O";
oLocations = locationData.trajectories(oIndices,:,:);
oColor = [77 240 214]/255; % cyan
oColor = [240 77 77]/255; % red
oSymbole = 'o';
oLineWidth = 6;

cIndices = atomNameArray == "C";
cLocations = locationData.trajectories(cIndices,:,:);
cColor = [240 77 77]/255; % red
cColor = [77 240 214]/255; % cyan
cSymbole = 'o';
cLineWidth = 6;

nIndices = atomNameArray == "S";
nLocations = locationData.trajectories(nIndices,:,:);
nColor = [44 0 107]/255; % violett
nSymbole = 'o';
nLineWidth = 1;

pIndices = atomNameArray == "N";
pLocations = locationData.trajectories(pIndices,:,:);
pColor = [0 0 0]/255; % black
pSymbole = 'o';
pLineWidth = 1;

sIndices = atomNameArray == "S";
sLocations = locationData.trajectories(sIndices,:,:);
sColor = [229 232 0]/255; % yellow
sSymbole = 'o';
sLineWidth = 1;


%% write Video
% close(writerObj);
writerObj = VideoWriter(sprintf('%s%s%s' ...
    ,resultsDir,datestr(now,"yyyymmdd_HHMM"),"_MDSimVideo"),'MPEG-4');
writerObj.FrameRate = 10;
open(writerObj);

initializeFigure('lineWidth',1.5,'legend',false,'posAndSize',[50 50 1920*1.5 1080*1.5]);
deltaViewAngle = 0.3;

for frameNr = 1:round(size(hLocations,3)*0.3)
    
    plotLocations(hLocations(:,:,frameNr),hColor,hSymbole,hLineWidth);
    hold on
    plotLocations(cLocations(:,:,frameNr),cColor,cSymbole,cLineWidth);
    plotLocations(oLocations(:,:,frameNr),oColor,oSymbole,oLineWidth);
    plotLocations(pLocations(:,:,frameNr),pColor,pSymbole,pLineWidth);
    plotLocations(sLocations(:,:,frameNr),sColor,sSymbole,sLineWidth);
    hold off
    axis([ -7 7 -5 5 -5 5]);
    view([-20+deltaViewAngle*frameNr 30]);
    xlabel("X [nm]");
    ylabel("Y [nm]");
    zlabel("Z [nm]");
    title(sprintf("Simulation duration: %.3f ns, deltaT: %.3f ps" ...
        ,simDurInNs,deltaTInPs));
    legend("H","C","O","P","S",'Location','bestoutside');
    
    writeVideo(writerObj,getframe(gcf));
end

close(writerObj);

%% convert to .mp4
% pathToVideo = writerObj.Path + "\" + writerObj.Filename;
% pathToVideoMP4 = regexprep(pathToVideo,'\.avi','.mp4'); 
% 
% [~,~] = system(sprintf('ffmpeg.exe -i %s -y -an -c:v libx264 -crf 0 -preset slow %s',pathToVideo,pathToVideoMP4)); % for this to work, you should have installed ffmpeg and have it available on PATH

% 
% %% test
% Z = peaks;
% surf(Z); 
% axis tight manual 
% set(gca,'nextplot','replacechildren'); 
% 
% v = VideoWriter('peaks','Uncompressed AVI');
% v.FrameRate = 5;
% open(v);
% 
% for k = 1:20 
%    surf(sin(2*pi*k/20)*Z,Z)
%    frame = getframe(gcf);
%    writeVideo(v,frame);
% end
% 
% close(v);

%% fumction

function plotLocations(locationsAtFrameX,atomColor,symbole,lineWidth)
plot3(locationsAtFrameX(:,1),locationsAtFrameX(:,2) ...
    ,locationsAtFrameX(:,3),symbole,'Color',atomColor ...
    ,'LineWidth',lineWidth);

end

function moleculeIndex = getMoleculeIndexFromMoleculeName(moleculeName)
indexAsCell = textscan(moleculeName,'%7d');
moleculeIndex = indexAsCell{1};

end








