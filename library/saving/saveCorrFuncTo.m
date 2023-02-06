function saveCorrFuncTo(corrFuncDirectory,corrFuncForDifferentNN,atomIndex ...
    ,NNcases,simDurInNs,deltaTInPs,theta,phi)

corrFuncOrder = strrep(inputname(2),"corrFunc","");
for numberOfNNCases = 1:size(corrFuncForDifferentNN,1)
    corrFunc = corrFuncForDifferentNN(numberOfNNCases,:); %#ok<NASGU>
    savingName = sprintf("%s%s_atomIndex%i_NN%i_simDur%sns_dT" ...
        + "%sps_theta%s_phi%s.mat",corrFuncDirectory,corrFuncOrder ...
        ,atomIndex,NNcases(numberOfNNCases) ...
        ,getDecimalNumberAsSeparatedWithUnderScore(simDurInNs) ...
        ,getDecimalNumberAsSeparatedWithUnderScore(deltaTInPs) ...
        ,getDecimalNumberAsSeparatedWithUnderScore(theta) ...
        ,getDecimalNumberAsSeparatedWithUnderScore(phi));
    if exist(savingName,'file')
        error("Correlation function already exists.");
    end
    save(savingName,'corrFunc','-v7.3')
end

if ~isfile(savingName)
    error("%s was not saved",savingName);
end

end