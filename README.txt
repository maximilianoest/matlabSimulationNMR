This is the README file to the existing Simulation.

Following information is important:
- the given mainMagnetigField value in the configuration file is only
    important to calculate a rough estimate of R1 and other constants. For
    better results, please use the calculated and summed correlation 
    functions to determine valid results.
- all set up parameters can be adjusted in the configMain.txt located at
    txtFiles directory. The path is hardcoded in simulationMain.m and is 
    maybe used somewhere else in the code.
- other scripts can be added to the scripts directory to simulate other 
    parameters. This script name should then be given in the configMain.txt
- nearestNeigbours is used if one case for NN should be used and 
    nearestNeighbourCases is used if multiple NN cases should be considered.
    But be careful: the script file should consider the right parameter
- in some files a RESULTS directory is used. This directory will not be 
    shipped with a github clone but can be found on the server and must be
    downloaded separately (too high data amount)
- values in constants.txt are important for the calculations. The path to
    the constants.txt file is hardcoded in simulation main and maybe is 
    used somewhere else in the code. 
- the trajectories are saved in single precision. Therefore, all data is 
    calculated in single precision, too. This is validated for the 
    calculation of the correlation function. There, the difference in 
    between double and single precsion is much smaller the the values of 
    the correlation function and its changed noise (calculated for 25 
    atoms and shown in saved under validateMethodsForComplexCorrFunc.m)
- for the calculation of the spectral density it is switched back to 
    double precision because the numbers calculated in the explicit 
    FT have some small and some large orders.


