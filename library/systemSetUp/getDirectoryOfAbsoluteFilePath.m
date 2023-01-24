function directory = getDirectoryOfAbsoluteFilePath(fileDirString)

splittedDirectoryPath = strsplit(fileDirString,filesep);
directory = strrep(fileDirString,splittedDirectoryPath{end},"");

end
