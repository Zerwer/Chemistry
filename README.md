# Chemistry
Neural networks for predicting chemical properties

# Octanol-Water Partition Coefficient 
Data: 24652 compounds scraped from PubChem

Regression predicts LogP from compounds Atom Pair fingerprint

See models/octanol_water_partition_coefficient/regression.py

# Water Solubility
Data: 9982 compounds from AqSolDB

Two regression models predict LogS, one from the Atom Pair fingerprint and the other from the predicted LogP. Both predictions are combined in another neural network to predict with higher accuracy. 

LogP see models/solubility_logP/regression.py  
Atom Pair see models/water_solubility/regression.py  
Combined see models/combined_water_solubility/regression.py  

# Boiling Point
Data: 2694 compounds scraped from PubChem and sorted

Classification model to determine best fingerprint for this application, predicts range of boiling point for a compound.

See models/boiling_point/simple_ia.py
