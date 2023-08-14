function saveFigureAsTiffTo(savingDirectory,whichLipid,simulationDate ...
    ,figureDescription,varargin)

if isempty(varargin)
   backUpOldFigure = true;
else
    backUpOldFigure = varargin{1};
end

if ~strcmp(savingDirectory(end),filesep)
    savingDirectory = sprintf('%s%s',savingDirectory,filesep);
end
if ~exist(savingDirectory,'dir')
    mkdir(savingDirectory)
end
fileName = sprintf('%s_%s_%s',whichLipid,simulationDate,figureDescription);
fileSavingPath = sprintf('%s%s.tiff',savingDirectory,fileName);

if exist(sprintf('%s%s',fileSavingPath),'file') && backUpOldFigure
    oldDirectory = sprintf('%sOLD%s',savingDirectory,filesep);
    if ~exist(oldDirectory,'dir')
        mkdir(oldDirectory);
    end
    movefile(fileSavingPath,sprintf('%s%s%s.png',oldDirectory ...
        ,datestr(now,'yyyymmdd_HHMM_'),fileName),'f');
end

print(gcf,fileSavingPath,'-dtiff','-r800');

end