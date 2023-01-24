function saveCorrFuncTo(corrFuncDirectory,corrFunc,atomIndex ...
    ,NNcases,simDurInNs,deltaTInPs)

orderOfCorrFunc = strrep(inputname(2),"corrFunc","");
for numberOfNNCases = 1:size(corrFunc,1)
    savingName = sprintf("%s%s_atomIndex%i_NN%i_simDur%sns_dT" ...
        + "%sps.mat",corrFuncDirectory,orderOfCorrFunc,atomIndex ...
        ,NNcases(numberOfNNCases) ...
        ,replace(num2str(simDurInNs),".","_") ...
        ,replace(num2str(deltaTInPs ),".","_"));
    save(savingName,'corrFunc','-v7.3')
end

if ~isfile(savingName)
    error("%s was not saved",savingName);
end

end