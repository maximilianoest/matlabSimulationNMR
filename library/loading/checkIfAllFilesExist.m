function checkIfAllFilesExist(filePaths)

if isempty(filePaths)
    error("No file paths were given.");
end

for filePath = filePaths
    if ~exist(filePath,'file')
        error("File: %s does not exist.",filePath);
    end
end
end