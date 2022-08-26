function createNewLogFileAndRenameOldOne(logFilePath)

if isfile(logFilePath)
    movefile(logFilePath,sprintf('%s_OLD%s.txt' ...
        ,strrep(logFilePath,'.txt',''),datestr(now,'HHMMSS')),'f');
% else
%     fileId = fopen(logFilePath);
%     fclose(fileId);
end

end
