clc; clear all; close all
addpath(genpath("..\molecularDynamicSimulationAnalysis"));
addpath(genpath("..\distributionAnalysis"));

run("tissueR1BasedOnLiterature");
run("determineR1OfMyelinBasedOnCompartments");
run("determineSurfaceWaterAndHistologR1InMyelin");
run("determineProteinR1WithRealisticFEModel");
run("fittingCompartmentR1ToPowerFunction");
run("NotFastExchangingModel");
