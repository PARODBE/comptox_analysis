[![Python](https://img.shields.io/pypi/pyversions/torchquad)](https://img.shields.io/pypi/pyversions/torchquad)
[![License](https://img.shields.io/badge/license-GPLv3-blue)](https://img.shields.io/badge/license-GPLv3-blue)

### 
``highlighting-atoms`` is a python package in order to highglight the fingerprints that the user selects, which can be useful for example if you obtain the importance of the variables through a Machine Learning model.

##### Install highlighting-atoms from PyPI
```bash
pip install highlighting-atoms
```

##### The following functions are available after installation:

```python
# Import library
from highlighting_atoms.highlighting_atoms import highlighter

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
