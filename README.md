#########################################
Prerequisites:
The code requires the following libraries
matplotlib.pyplot 
numpy 
scipy.integrate 
from scipy.interpolate import interp1d
time
##########################################
plotly.graph_objects - for plotting the map with the stages downranges
To be able to plot the maps in the code you need a mapbox token which can be obtained here after registration: https://docs.mapbox.com/help/glossary/access-token/
The maps are not essential for the results from the code so you can just comment them. You need to comment just the function plotmap().


HOW TO WORK WITH THE CODE

input.py is the input file for Saturn_Model.py
inputElectron.py is the input file for Eletron_Model.py
inputNovaRocket.py is the input file for NovaRocket.py

To run any of these models you need to set up the parameters in the input file and then run the model file that you want to.
