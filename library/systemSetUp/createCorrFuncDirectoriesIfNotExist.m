function corrFuncDirs = createCorrFuncDirectoriesIfNotExist( ...
    resultsDirectory,gromacsSimulationDate)

if ~ endsWith(resultsDirectory,filesep)
    resultsDirectory = resultsDirectory + filesep;
end

corrFuncDirs = gromacsSimulationDate + "_corrFunc" ...
        + ["ZerothOrder" "FirstOrder"  "SecondOrder"] + filesep;

for directory = corrFuncDirs
    if ~isfolder(resultsDirectory+ directory)
        mkdir(resultsDirectory + directory)
    end
end

end