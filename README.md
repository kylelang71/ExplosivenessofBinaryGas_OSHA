# ExplosivenessofBinaryGas_OSHA
This program was written to analyse simple industrially available chemical mixtures used as ingredients in consumer products.  It determines and reports whether the vapor phase in equilibrium with a given binary fluid is potentially explosive as defined by OSHA at ambient temperature and pressure.  The Peng-Robinson equation of state was used to model the compostion of the vapor.  This equation of state is not a useful predictor of the properties of all chemical species, familiarize yourself with the model and its shortcomings before use. 
main.py 
  Main script.  Given the composition of a binary liquid, physical properties of each component are pulled into a dataframe from larger datasets. It has logic to operate on the output of PR__EOS.py and provide information regarding the explosive potential of the vapor. 
PR_EOS.py
  Function script.  This script is called and passed the relevant values in the dataframe to predict vapor composition at ambient temperature and pressure. The composition of the vapor is passed back to main.py.
Explosive_Limits.csv
  OSHA sourced dataset of Upper and Lower Explosive limits of common pure chemicals.  Future update - ENSURE THAT THE COMPOUNDS OF INTEREST ARE INCLUDED IN THIS DATASET TO AVOID ERRORS.  
Physical_Properties.csv
  Pure chemical physical property data from 'Introductory Chemical Engineering and Thermodynamics' Second Edition by Prentice Hall / Elliot & Lira.  Future update - ENSURE THAT THE COMPOUNDS OF INTEREST ARE INCLUDED IN THIS DATASET TO AVOID ERRORS.  
Glass_Cleaner.csv
  An ammonium - water mixture, used as a baseline efficacy test for glass cleaners.  This is a placeholder used as an example of the program's capabilities. Alter this for further investigations.  
