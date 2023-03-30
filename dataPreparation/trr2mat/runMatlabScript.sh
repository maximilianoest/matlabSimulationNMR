#!/bin/bash

matlab -nosplash -nodisplay -nodesktop -r 'try; createMatFileFromTerminateTrrFile; catch e; disp(e.message); disp(e.stack); save(sprintf("%s_Error",datestr(now,"yyyymmdd")),"e");  end; quit;' > Bilayer_distributionAnalysis.log