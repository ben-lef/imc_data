# Supplementary data for the articles:
## 1) "Quantification of intermetallic compound layer morphology in aluminum-steel arc welding"
## 2) "Simulation of the intermetallic compound layer formation in aluminum-steel arc welding"

The Python file `code.py` is furnished to process the dataset taining all the information about the IMC reported in the aforementioned articles. 
The dataset `imc_data.csv` is a large file and could not be directly stored in GitHub. It can be downloaded from the following link: https://leflon.fr/benjamin/imc_data.csv

To immediately put the code and database together in application, make sure to **store both** `code.py` **and** `imc_data.csv` **in a same directory.**  
As an example, you could try the following call in the environment of `code.py`:

 `full_analysis(db, 'i140v140', 1, 'T30')`
 
As is, this will produce various plots describing the morphology of the layer for the i140v140 condition, on its first experimental attempt, on a transverse cross-section extracted at z = 30 mm. 
More importantly, it will return a dictionnary containing the ensemble of the local, global, micro- and macro-descriptors as defined in the article, for said sample.
To visualize the profile of any local descriptor along the baseline of the cross-section, it should be plotted against the list of positions which is stored behind the keyword "pos" in the dictionnary.

Note: the dataset used in article 2) is an extension from the dataset used in article 1). The former contains additional experimental data, with a new welding condition (i150v230) and data obtained when using a copper backing plate rather than a steel one. Data involving the use of a copper backing plate have the prefix "cu_" in the label column. For instance, the label "i150v230" identifies welds obtained with a steel backing plate, while "cu_i150v230" identifies welds obtained with similar parameters but a copper backing plate.
