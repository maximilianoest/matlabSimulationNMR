classdef generalSimulation
    %GENERALSIMULATION Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        fileDirectory;
        fileName;
        runOnServer;
        trajectoriesX;
        trjaectoriesY;
        trajectoriesZ;
        gromacsDatasetConfiguration;
        dependencies;
        
    end
    
   
    
    methods
        function this = generalSimulation(fileDirectory,fileName)
            determineRunOnServer(this);
            this.fileDirectory = fileDirectory;
            this.fileName = fileName;
            this.FilePath = sprintf('%s%s%s',fileDirectory,filesep, ...
                fileName)
            setUpDependencies(this);
        end
        
        
    
        function loadGromacsDataFromMatFile(this)
            matObject = matfile(path2File);
            simConfig = matObject.configuration;
            
            load(path2File,'trajectories');
            
            trajectoryX = squeeze(single(trajectories(:,1,:)));
            trajectoryY = squeeze(single(trajectories(:,2,:)));
            trajectoryZ = squeeze(single(trajectories(:,3,:)));
            
        end
        
        function getGromacsSimulationConfiguration(this)
           
            
        end
        
        function logGromacsSimulationConfiguration(this)
        
            
        end
        
        function saveResults(this)
            createDataSavingObject();
            
        end
    end
    
    methods (Access = protected)
        
        function runOnServer = determineRunOnServer()
            plattform = computer;
            switch plattform
                case {'PCWIN64'}
                    runOnServer = 0;
                case {'GLNXA64'}
                    runOnServer = 1;
                otherwise 
                    error('generalSimulation:determineRunOnServer', ...
                        'plattform: %s is not part of plattform list', ...
                        plattform);
            end
        end
        
%         function setUpDependencies(this)
%             
%             lipidName = getLipidNameFromFileName(this.fileName);
%             constituent = getConstituentFromFileName(fileName);
%             combinedName = [lipidName constituent];
%             
%             runOnServer = configuration.runOnServer;
%             if runOnServer
%                 path2BaseConfiguration = configuration.path2BaseConfigurationOnServer;
%             else
%                 path2BaseConfiguration = ...
%                     configuration.path2BaseConfigurationOnLocalMachine;
%             end
%             
%             baseConfiguration = readConfigurationFile(path2BaseConfiguration);
%             startingDate = datestr(date,'yyyymmdd');
%             
%             if runOnServer
%                 path2Data = baseConfiguration.path2DataOnServer;
%                 path2ConstantsFile = baseConfiguration.path2ConstantsFileOnServer;
%                 mainPath2Results = baseConfiguration.path2ResultsOnServer;
%                 additionalFolder = '';
%             else
%                 path2Data = baseConfiguration.path2DataOnLocalMachine;
%                 path2ConstantsFile = baseConfiguration.path2ConstantsFileOnLocalMachine;
%                 mainPath2Results = baseConfiguration.path2ResultsOnLocalMachine;
%                 additionalFolder = [configuration.kindOfResults '\'];
%             end
%             
%             path2Results = sprintf('%s%s',mainPath2Results,additionalFolder);
%             if ~isfolder(path2Results)
%                 mkdir(path2Results);
%             end
%             
%             path2Save = sprintf('%s%s_Results_%s_%s.mat',path2Results,startingDate ...
%                 ,configuration.kindOfResults,combinedName);
%             path2LogFile = sprintf('%s%s_LogFile_%s_%s.txt',path2Results ...
%                 ,startingDate,configuration.kindOfResults,combinedName);
%             
%             createNewLogFileAndRenameOldOne(path2LogFile);
%             
%         end
    end
    
     methods(Abstract)
        loadConfiguration(this);
        performCalculations(this);
        createDataSavingObject(this);
    end
end

