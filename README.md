[![Python](https://img.shields.io/pypi/pyversions/torchquad)](https://img.shields.io/pypi/pyversions/torchquad)
[![License](https://img.shields.io/badge/license-GPLv3-blue)](https://img.shields.io/badge/license-GPLv3-blue)

### 
In January 2024, I defended my doctoral thesis in computational toxicology, and since then, I have been contemplating how to contribute to the community with code that can simplify calculations needed according to the specific problem. Over these years, among many other things, I have worked on topics related to interpretability, visualization of chemical space, and combining the outputs of different models. In all cases, the problem was common: I often needed a lot of code for tasks that were part of my daily routine. That's why I decided to create ``comptox_analysis``. In the case of visualizing the chemical space, if molecular descriptors are selected, a PCA will be represented and in the case of fingerprints, a t-SNE. Fingerprint highlighting is compatible with both morgan and rdkit fingerprints.


##### Install highlighting-atoms from PyPI
```bash
pip install comptox_analysis
```

##### The following functions are available after installation:

```python
# Import library
from highlighting_atoms.visualizing import chemical_space_plotter
from highlighting_atoms.ML_stacking import Stacking
from highlighting_atoms.highlighting import highlighter
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC
from sklearn.linear_model import LogisticRegression
from sklearn.neighbors import KNeighborsClassifier

chemical_space_plotter('cox2_train.sdf').visualizer(type='morgan',ID_column='pchembl_value')

models_params_data = {
    'RandomForest': (RandomForestClassifier(random_state=42), 
                     {'n_estimators': [50, 100], 'max_depth': [None, 10]}, 
                     (X2_train, y2_train)),
    'SVM': (SVC(probability=True, random_state=42), 
            {'C': [0.1, 1, 10], 'kernel': ['linear', 'rbf']}, 
            (X3_train, y3_train)),
    'LogisticRegression': (LogisticRegression(random_state=42), 
                           {'C': [0.1, 1, 10]}, 
                           (X4_train, y4_train)),
    'KNN': (KNeighborsClassifier(), 
            {'n_neighbors': [3, 5, 7]}, 
            (X5_train, y5_train))
}

# # Inicializar el modelo de Stacking
stacking_model = Stacking(models_params_data=models_params_data).fit(X_test, y_test)

params = {
    'logical_rule': ['OR', 'AND', 'Majority']
}

# Realizar GridSearchCV en el modelo de Stacking
grid_search_stacking = GridSearchCV(stacking_model, params, cv=5, scoring='accuracy')
grid_search_stacking.fit(X_test, y_test)

#for one compound
mols=['CC(C)(C)NC(=O)[C@@H]1C[C@@H]2CCCC[C@@H]2CN1C[C@@H](O)[C@H](Cc1ccccc1)NC(=O)[C@H](CC(N)=O)NC(=O)c1ccc2ccccc2n1',
 'C/C=C/C[C@@H](C)[C@@H](O)[C@H]1C(=O)N[C@@H](CC)C(=O)N(C)CC(=O)N(C)[C@@H](CC(C)C)C(=O)N[C@@H](C(C)C)C(=O)N(C)[C@@H](CC(C)C)C(=O)N[C@@H](C)C(=O)N[C@H](C)C(=O)N(C)[C@@H](CC(C)C)C(=O)N(C)[C@@H](CC(C)C)C(=O)N(C)[C@@H](C(C)C)C(=O)N1C',
 'CC1(C)S[C@@H]2[C@H](NC(=O)[C@H](N)c3ccc(O)cc3)C(=O)N2[C@H]1C(=O)O',
 'CC(C)(C)S(=O)(=O)C[C@@H](Cc1ccccc1)C(=O)N[C@@H](Cc1c[nH]cn1)C(=O)N[C@@H](CC1CCCCC1)[C@@H](O)[C@@H](O)C1CC1',
 'O=C(N[C@H](CO)[C@H](O)c1ccc([N+](=O)[O-])cc1)C(Cl)Cl']

highlighter_instance = highlighter([mols[0]]) #this is for one compound

highlighter_instance.render_image(number=0,indexes=False)

highlighter_instance.fragmentation(n=10,number=0)[0] # [0] is to visualize the fragments and [1] to obtain the fragment dict.


highlighter_instance.highlighting(type='morgan',fingerprint_numbers=[[707]]) #this is for one compound. You can specify the fp type and its hyperparameters

#for multiple compounds

highlighter_instance = highlighter(mols)

highlighter_instance.highlighting(fingerprint_numbers=[[3,779],[3,779],[3,779],[3,779],[3,779]])
