# Java_Analysis_Suite
Java Analysis Suite for CLAS12 Data

The Analysis pipline code is in the Phd_Clean folder. The code clean the raw data acquired during Fall2018 and produce a new final file of about size 1/100 of original data including some reconstructed information regarding the interaction such as momenta, angle, transverse momentum , pseudo-rapidity, hadron angle in the COM frame.

The other folders provide compiled services such as creating histograms, creating visualization tools, providing function minimizations for fitting and regressions, etc..

Minimization use the FreeHEP Minuit package : http://java.freehep.org/freehep-jminuit/apidocs/index.html
Particle class, Filter classes and Histogram classes are from the standard JNP and JROOT classes we developed at Jlab (the source can be found in the other repositories on my GitHub).


1) create a jar file using the given Manifest.
2) Use the auger script file present on the JLAB ifarm at : /work/clas12/gangel/phd/ 
3) Submit the script to the farm linking the jar.
4) The code will run on all the data acquired in Fall 2018, and create a new Hipo file with reconstructed information.

