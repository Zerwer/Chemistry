# Chemistry
Neural networks for predicting chemical properties and UI for graphing predict properties of datasets.
Developing generative adversarial network that learns to create molecules based on configured predicted properties 
and other constants such as molecular weight.

## Conda Setup

```
conda create -n Chemistry
conda activate Chemistry
conda install -c rdkit rdkit
conda install scikit-learn
conda install matplotlib
```

# Interactive UI
Requires models to be compiled and saved.
Molecules can be manually entered as SMILES strings or a file can be loaded.

## Main Window 
- Table containing all loaded molecules
- Search bar to search through molecules
- Property display window to show predicted properties of a selected molecule

## Graph Generation
Select graph type from View>Graph tab, a window to configure the graph then appears.

# Property Prediction

## Octanol-Water Partition Coefficient (logP)
Data: 24652 compounds scraped from PubChem

Regression predicts LogP from compounds Atom Pair fingerprint

See models/octanol_water_partition_coefficient/regression.py

## Water Solubility
Data: 9982 compounds from AqSolDB

Two regression models predict LogS, one from the Atom Pair fingerprint and the other from the predicted LogP. Both predictions are combined in another neural network to predict with higher accuracy. 

LogP see models/solubility_logP/regression.py  
Atom Pair see models/water_solubility/regression.py  
Combined see models/combined_water_solubility/regression.py  

## Melting Point
Data: 7952 compounds from PubChem

Regression model predicts melting point from logP model and water solubility.

See models/melting_point/reverse_gse_regression.py

## Boiling Point
Data: 2694 compounds scraped from PubChem and sorted

Classification model to determine best fingerprint for this application, predicts range of boiling point for a compound.

See models/boiling_point/simple_ia.py

## Data Collection
To run the models you will need data to train the models.

For data that collects from PubChem go to: https://pubchem.ncbi.nlm.nih.gov/classification/  
1. Under select classification select PubChem>PubChem Compound TOC.
2. In the tree select Chemical Properties>Experimental Properties.
3. Click on the number next to property and open new page.
4. On the right click structure download.
5. Select SMILES format and leave the rest and download.
6. Use this file as input list for data_collection

### LogP/Octanol Water Partition Coefficient
Run data_collection/octanol_water_partition_coefficient.py and select the path to the PubChem list of compounds(See above)

### Solubility
Download https://www.amdlab.nl/database/AqSolDB/

Run data_collection/aqsoldb_transform.py on the csv to change the format.

### Melting Point
Run data_collection/melting_point.py and select the path to the PubChem list of compounds(See above)

