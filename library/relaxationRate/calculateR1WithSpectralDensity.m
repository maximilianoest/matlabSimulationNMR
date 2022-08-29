function [r1WithPerturbationTheory] = calculateR1WithSpectralDensity(w1 ...
    ,w2,DD)
% This function calculates the R1 Rate for given spectral densities

r1WithPerturbationTheory = DD*3/2*(abs(real(mean(w1))) ...
    +abs(real(mean(w2))));
    
end
