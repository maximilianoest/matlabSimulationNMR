classdef lipid < handle
    properties
        lipidName
        atomComposition
        lipidWeight
        lipidClass
        relativeFrequencyInMyelinModel
        numberOfMoleculesInModel = 1
        frequencyOfLipidWithSameClass
    end
    
    methods
       function obj = lipid(lipidName,atomNameCollectionLength)
           obj.lipidName = lipidName;
           obj.atomComposition = zeros(1,atomNameCollectionLength);
       end
       
       function increaseAtomCompositionCounter(obj,indexToIncrease)
          obj.atomComposition(indexToIncrease) = ... 
              obj.atomComposition(indexToIncrease) + 1;
       end
       
       function calculateLipidWeight(obj,atomWeightCollection)
           obj.lipidWeight = obj.atomComposition * atomWeightCollection';
       end
       
       function determineLipidClass(obj,lipidComposition)
           try
               indices = obj.getIndicesForClassDetermination(obj.lipidName, ...
                   lipidComposition);
           catch exception
               if exception.identifier == "LIPID:tooManyIndices" && ...
                       obj.lipidName == "PSPI"
                   indices = obj.getIndicesForClassDetermination( ...
                      extractAfter(obj.lipidName,"PS"),lipidComposition);
               else 
                   rethrow(exception)
               end
           end
           obj.lipidClass = lipidComposition(logical(indices));
       end
           
       function indices = getIndicesForClassDetermination(~, ...
               nameOfLipid,lipidComposition)
           indices = zeros(1,length(lipidComposition));
           
           for lipidClassNr = 1:length(lipidComposition)
               if contains(nameOfLipid,lipidComposition(lipidClassNr))
                   indices(lipidClassNr) = 1;
               end
           end
           
           if sum(indices) < 1
               error("LIPID:tooLessIndices" ...
                   ,"No appropriate lipid class was found. (%s)" ...
                   ,nameOfLipid);
           elseif sum(indices) > 1
               error("LIPID:tooManyIndices" ...
                   ,"More than one lipid class was found. (%s)" ...
                   ,nameOfLipid);
           end
       end
       
       function determineRelativeFrequencyInMyelinModel(obj ...
        ,lipidComposition,averageLipidDistribution)
            obj.relativeFrequencyInMyelinModel = averageLipidDistribution( ...
                lipidComposition == obj.lipidClass);
       end
    end
end
    

