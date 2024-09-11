[![Python](https://img.shields.io/pypi/pyversions/torchquad)](https://img.shields.io/pypi/pyversions/torchquad)
[![License](https://img.shields.io/badge/license-GPLv3-blue)](https://img.shields.io/badge/license-GPLv3-blue)

### 
In January 2024, I defended my doctoral thesis in computational toxicology, and since then, I have been contemplating how to contribute to the community with code that can simplify calculations needed according to the specific problem. Over these years, among many other things, I have worked on topics related to interpretability, visualization of chemical space, and combining the outputs of different models. In all cases, the problem was common: I often needed a lot of code for tasks that were part of my daily routine. That's why I decided to create ``comptox_analysis``.


# Installation

## Dependencies

comptox_analysis requires:

- rdkit=2023.09.6
- bokeh=3.1.1
- scikit-learn
- matplotlib
- numpy
- pandas

_ _ _

## User installation

The easiest way to install comptox_analysis is using ```pip```:

```bash
pip install comptox_analysis
```

_ _ _

# Example of Usage for Interpretability using highlighter

Imagine you obtain the feature importance, which in this case are fingerprints, using a ML or DL model, and you want to highlight the most important fingerprints for each molecule. To carry out this task, we first load the necessary libraries and the dataset, build a Random Forest model as an example, and then extract the feature importance. 

```python
# Import library
from comptox_analysis.highlighting.highlighting_atoms import highlighter
from rdkit import Chem
from rdkit.Chem import PandasTools
from rdkit.Chem.AllChem import GetMorganFingerprintAsBitVect
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestRegressor

df=PandasTools.LoadSDF('caco2.sdf',smilesName='Smiles')

def fingerprints_inputs2(dataframe):
        X=np.array([GetMorganFingerprintAsBitVect(mol,radius=2,nBits=2048,useFeatures=True) for mol in [Chem.MolFromSmiles(m) for m in list(dataframe.Smiles)]])
        y=dataframe.activity.astype('float')
        return X,y

X,y=fingerprints_inputs2(df)

model=RandomForestRegressor(random_state=46).fit(X,y)

columns=[f'fp_{i}' for i in range(2048)]

imp=pd.Series(data=model.feature_importances_,index=columns).sort_values(ascending=False)
imp[:10].plot.bar();

```

<p align="center">
  <img src="https://raw.githubusercontent.com/phi-grib/comptox_analysis/main/images/importances.PNG" alt="Cover Page">
</p>

We are going to select compounds with two of the ten most important variables (in this case, we obtain a list with three compounds). Now, we render the image of the first compound using ```render_image``` function:

```python

df_fp=pd.DataFrame(X,columns=columns)

indexes=df_fp[(df_fp.fp_576==1)&(df_fp.fp_779==1)].index.values

mols=df.loc[indexes,'Smiles'].tolist()

highlighter_instance=highlighter(mols)

highlighter_instance.render_image(number=0,indexes=False)

```

<p align="center">
  <img src="https://raw.githubusercontent.com/phi-grib/comptox_analysis/main/images/molecule_rend..PNG" alt="Cover Page">
</p>


Also, you can visualize molecular fragments using ```fragmentation``` function:

```python

highlighter_instance.fragmentation(n=26,number=0)[0]

```

<p align="center">
  <img src="https://raw.githubusercontent.com/phi-grib/comptox_analysis/main/images/frag..PNG" alt="Cover Page">
</p>

And now, using ```highlighting``` function, you will be able to observe the fragments that you want to highlight (the package supports both, rdkit and morgan fingerprints):

```python

highlighter_instance = highlighter([mols[0]])

highlighter_instance.highlighting(type='morgan',fingerprint_numbers=[[1300]])

```

<p align="center">
  <img src="https://raw.githubusercontent.com/phi-grib/comptox_analysis/main/images/one_compound_highl..PNG" alt="Cover Page">
</p>

Also, it is possible to highlight multiple compounds (in this case we are highlighting the selected fingerprints as an example9:

```python

#for multiple compounds

highlighter_instance = highlighter(mols)

highlighter_instance.highlighting(fingerprint_numbers=[[576,779],[576,779],[576,779]])

```

<p align="center">
  <img src="https://raw.githubusercontent.com/phi-grib/comptox_analysis/main/images/multiple_highl..PNG" alt="Cover Page">
</p>
