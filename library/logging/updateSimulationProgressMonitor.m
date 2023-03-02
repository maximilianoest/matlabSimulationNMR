function updateSimulationProgressMonitor(corrFuncZerothOrder ...
    ,corrFuncFirstOrder,corrFuncSecondOrder,avgSpecDens1,avgSpecDens2 ...
    ,r1Estimation,dTInSec,avgRegionSpecDens,resultsDirectory ...
    ,whichLipid,matlabSimulationDate)

fig = initializeFigure();
initializeSubplot(fig,2,2,1);
timeAxis = 0 : dTInSec : (length(corrFuncZerothOrder) -1 ) * dTInSec;
plotCorrFuncs(corrFuncZerothOrder,corrFuncFirstOrder ...
    ,corrFuncSecondOrder,timeAxis);
title("Correlation functions");
legend("real(0th)","imag(0th)","real(1st)","imag(1st)" ...
    ,"real(2nd)","imag(2nd)");
xlabel("lag time [s]");

initializeSubplot(fig,2,2,2);
title("Real(spectral density)");
specDensTimeAxis = timeAxis(round(end*avgRegionSpecDens(1) ...
    :round(end*avgRegionSpecDens(2))));
plot(specDensTimeAxis,real(avgSpecDens1));
plot(specDensTimeAxis,real(avgSpecDens2));
legend("1st", "2nd");
xlabel("integration limit [s]");

initializeSubplot(fig,2,2,3);
title("r1Estimation");
calculatedAtomsAxis = 1:length(r1Estimation);
plot(calculatedAtomsAxis,r1Estimation)
axis([ 0 inf 0 inf])
leg = legend();
leg.Visible = "off";

ax = initializeSubplot(fig,2,2,4);
fourthSubplotText = sprintf("Current R1 estimate: %.4f",r1Estimation(end));
text(0.05,0.5,fourthSubplotText,'interpreter','latex');
set(ax,'visible','off');
lgd = legend();
lgd.Visible = 'Off';
simDuration = length(corrFuncFirstOrder)*dTInSec/1e-9;
dTInPs = dTInSec/1e-12;
specification = sprintf("simulationProgress_simDur%sns_dT%sps" ...
    ,getDecimalNumberAsSeparatedWithUnderScore(simDuration) ...
    ,getDecimalNumberAsSeparatedWithUnderScore(dTInPs));
saveFigureTo(resultsDirectory,whichLipid,matlabSimulationDate ...
    ,specification,false)

close(fig);

end

function plotCorrFuncs(corrFunc0,corrFunc1,corrFunc2,timeAxis)

plotRealAndImagPartOfCorrFunc(corrFunc0,timeAxis);
plotRealAndImagPartOfCorrFunc(corrFunc1,timeAxis);
plotRealAndImagPartOfCorrFunc(corrFunc2,timeAxis);

end

function plotRealAndImagPartOfCorrFunc(corrFunc,timeAxis)

plot(timeAxis,real(corrFunc));
plot(timeAxis,imag(corrFunc));


end