function saveCorrFuncTo(corrFuncDirectory,corrFuncForDifferentNN,atomIndex ...
    ,NNcases,simDurInNs,deltaTInPs,theta,phi)

orderOfCorrFunc = strrep(inputname(2),"corrFunc","");
for numberOfNNCases = 1:size(corrFuncForDifferentNN,1)
    corrFunc = corrFuncForDifferentNN(numberOfNNCases,:); %#ok<NASGU>
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