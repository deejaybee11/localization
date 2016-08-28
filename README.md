# Localization simulation

Simulates the expansion of a BEC inside a dumbell shaped potential with point scatterers along the channel.

Uses Normalized Gradient Flow (Also known as Imaginary Time Evolution) via the Time-Splitting Spectral Method (TSSP)http://www.aimsciences.org/journals/displayArticlesnew.jsp?paperID=8003 to calculate the ground state and in regular time to calculate the time evolution.

To run code on the NeSI cluster, you can update the files in your folder by running "$ bash update_cluster". Edit "submit.sl" to have the required number of cores, memory and time. Load the required modules (CFITSIO and intel compiler) by calling "$ ml CFITSIO" and "$ ml intel/2015a".
Submit the job with "$ sbatch submit.sl" 
Outputs are save periodically into the folder fits/date/.

Important parameters are found in "simulationdata.cpp" and are (supposed to be) located in similar sections and separated by comments. The number of data points used in the simulation are set when calling the SimulationData initializer function in the main. If you wish to profile the potential without building the simulation and looking at the dumbell each time you can use the python script included to generate an image as it uses the same parameters.
