function updateSimulationProgressForCrossAndAutoCorrelation( ...
    corrFuncZerothOrder,corrFuncFirstOrder,corrFuncSecondOrder ...
    ,avgSpecDensZero,avgSpecDensFirst,avgSpecDensSecond,dTInSec ...
    ,avgRegionSpecDens,varargin)



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
plot(specDensTimeAxis,real(avgSpecDensZero));
plot(specDensTimeAxis,real(avgSpecDensFirst));
plot(specDensTimeAxis,real(avgSpecDensSecond));
legend("0th", "1st", "2nd");
xlabel("integration limit [s]");

initializeSubplot(fig,2,2,3);
if isempty(varargin)
    title("No data for R1 are plotted.")
else
    title("Estimated R1");
    calculatedAtomsAxis = 1:length(varargin{1});
    lgdEntries = {};
    r1Estimates = "";
    for argInNr = 1:length(varargin)
        r1Name = inputname(nargin-length(varargin)+argInNr);
        lgdEntries{end+1} = r1Name; %#ok<AGROW>
        plot(calculatedAtomsAxis,varargin{argInNr});
        r1Estimates = sprintf("%s  %s: %.4f \n",r1Estimates,r1Name ...
            ,varargin{argInNr}(end));
    end
    axis([ 0 inf 0 inf])
    legend(lgdEntries);
    
    
    ax = initializeSubplot(fig,2,2,4);
    
    fourthSubplotText = sprintf("Current R1 estimates:\n%s" ...
        ,r1Estimates);
    text(0.05,0.5,fourthSubplotText,'interpreter','latex');
    set(ax,'visible','off');
    lgd = legend();
    lgd.Visible = 'Off';
end

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