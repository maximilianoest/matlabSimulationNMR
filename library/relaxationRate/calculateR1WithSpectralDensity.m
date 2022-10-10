function [r1WithPerturbationTheory] = calculateR1WithSpectralDensity(spectralDensityFirstOrder ...
    ,spectralDensitySecondOrder,dipoleDipoleConstant)
% This function calculates the R1 Rate for given spectral densities

r1WithPerturbationTheory = dipoleDipoleConstant*3/2*(abs(real(mean( ...
    spectralDensityFirstOrder))) + abs(real(mean( ...
    spectralDensitySecondOrder))));
    
end
