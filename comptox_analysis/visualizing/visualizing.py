import rdkit
from rdkit import Chem
from rdkit.Chem import PandasTools,Draw, AllChem, SDWriter
from rdkit.Chem import Descriptors
import pandas as pd
import numpy as np
from IPython.display import SVG
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem.PandasTools import ChangeMoleculeRendering
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from comptox_analysis.visualizing.predictable_tsne import PredictableTSNE
from sklearn.manifold import TSNE
from sklearn.neural_network import MLPRegressor
from bokeh.plotting import figure, show, output_notebook, ColumnDataSource
from bokeh.models import HoverTool
from bokeh.transform import factor_cmap
from bokeh.palettes import Spectral4
from bokeh.palettes import Category20_20
from bokeh.plotting import figure, output_file, save
from bokeh.palettes import RdBu3
from itertools import cycle
from bokeh.palettes import Category20
from typing import Tuple, List, Dict

class chemical_space_plotter:

    def __init__(self, path: str = None, test_path: str = None):
        '''
        This library allows drawing the chemical space for a specific dataset, providing structure chemical information.

        path corresponds to training set path and test_path to the test series.

        ----------------

        Example of use:
        chemical_space_plotter('train.sdf', 'test.sdf').visualizer(type='morgan', ID_column='pchembl_value', test=True)
        '''

        self.df = PandasTools.LoadSDF(path, smilesName='Smiles', includeFingerprints=True, molColName='Structure')
        self.df_test = PandasTools.LoadSDF(test_path, smilesName='Smiles', includeFingerprints=True, molColName='Structure') if test_path else None

    def _prepareMol(self, mol, kekulize):
        '''Prepare a molecule for visual representation by performing specific preprocessing steps.'''
        mc = Chem.Mol(mol.ToBinary())
        if kekulize:
            try:
                Chem.Kekulize(mc)
            except:
                mc = Chem.Mol(mol.ToBinary())
        if not mc.GetNumConformers():
            rdDepictor.Compute2DCoords(mc)
        return mc

    def moltosvg(self, mol, molSize=(450, 200), kekulize=True, drawer=None, **kwargs):
        '''Convert molecules into SVG format for visual representation.'''
        mc = self._prepareMol(mol, kekulize)
        if drawer is None:
            drawer = rdMolDraw2D.MolDraw2DSVG(molSize[0], molSize[1])
        drawer.DrawMolecule(mc, **kwargs)
        drawer.FinishDrawing()
        svg = drawer.GetDrawingText()
        return SVG(svg.replace('svg:', ''))

    def generate_colors(self, num_colors):
        '''Generate a set of distinct colors for different classes.'''
        palette = Category20[20]
        color_cycle = cycle(palette)
        colors = [next(color_cycle) for _ in range(num_colors)]
        return colors

    def visualizer(self, type='morgan', hyp_fp: Dict = None, activity_column='Activity', ID_column=None, test=False, Save=False, path='chemical_space_plot.html'):
        '''
        Visualize the chemical space of a dataset. If in type you select molecular_descriptor a PCA is built, however if you select morgan or rdkit a parametric TSNE is built training the TSNE output with a MLPRegressor. 
        hyp_fp allows to select the hyperparameters for the fingerprint generation.
        activity column corresponds to the name of the output variable. ID_column can be the molecule names or even pchembl_value. It is possible
        to save the HTML specifying Save option as True (if you want, specify the path). If test=True, a separate test set will also be visualized.
        '''
        # Load the training set
        mols_train = [Chem.MolFromSmiles(m) for m in self.df.Smiles]
        y_train = self.df[activity_column]

        if hyp_fp is None:
            if type == 'morgan':
                hyp_fp = {'radius': 2, 'nBits': 1024, 'useFeatures':True}  # Valores por defecto para Morgan fingerprints
            elif type == 'rdkit':
                hyp_fp = {'maxPath':2, 'fpSize': 2048}  # Valores por defecto para RDKit fingerprints

        if type == 'molecular_descriptor':
            df_molecular_descriptors_train = pd.DataFrame([Descriptors.CalcMolDescriptors(m) for m in mols_train])
            sc = StandardScaler().fit(df_molecular_descriptors_train)
            X_train_sc = sc.transform(df_molecular_descriptors_train)
            emb = PCA(n_components=2, random_state=46).fit(X_train_sc)
            emb_train = emb.transform(X_train_sc)

        elif type == 'morgan':
            X_fp_train = np.array([AllChem.GetMorganFingerprintAsBitVect(mol, **hyp_fp) for mol in mols_train])
            emb = PredictableTSNE(estimator=MLPRegressor(random_state=46), transformer=TSNE(random_state=46)).fit(X_fp_train, y_train)
            emb_train = emb.transform(X_fp_train)

        elif type == 'rdkit':
            X_fp_train = np.array([AllChem.RDKFingerprint(mol, **hyp_fp) for mol in mols_train])
            emb = PredictableTSNE(estimator=MLPRegressor(random_state=46), transformer=TSNE(random_state=46)).fit(X_fp_train, y_train)
            emb_train = emb.transform(X_fp_train)

        # Convert training embeddings to a DataFrame
        df_train_emb = pd.DataFrame(np.c_[emb_train, y_train], index=self.df.index.values, columns=['dim_1', 'dim_2', 'Activity'])

        # Test set processing (if test=True and a test dataframe is provided)
        if test and self.df_test is not None:
            mols_test = [Chem.MolFromSmiles(m) for m in self.df_test.Smiles]
            y_test = self.df_test[activity_column]

            if type == 'molecular_descriptor':
                df_molecular_descriptors_test = pd.DataFrame([Descriptors.CalcMolDescriptors(m) for m in mols_test])
                X_test_sc = sc.transform(df_molecular_descriptors_test)
                emb_test = emb.transform(X_test_sc)

            elif type == 'morgan':
                X_fp_test = np.array([AllChem.GetMorganFingerprintAsBitVect(mol, **hyp_fp) for mol in mols_test])
                emb_test = emb.transform(X_fp_test)

            elif type == 'rdkit':
                X_fp_test = np.array([AllChem.RDKFingerprint(mol, **hyp_fp) for mol in mols_test])
                emb_test = emb.transform(X_fp_test)

            # Convert test embeddings to a DataFrame
            df_test_emb = pd.DataFrame(np.c_[emb_test, y_test], index=self.df_test.index.values, columns=['dim_1', 'dim_2', 'Activity'])

            # Asignar un color único al conjunto de test (diferente de los colores de las clases)
            test_color = '#0000FF'  # Azul

        # Create graphical representation for training set
        unique_activities_train = y_train.unique()
        class_mols_train = {activity: [Chem.MolFromSmiles(smiles) for smiles in self.df[self.df[activity_column] == activity]['Smiles']] for activity in unique_activities_train}
        class_svgs_train = {activity: [self.moltosvg(mol).data for mol in mols] for activity, mols in class_mols_train.items()}
        
        sources_train = {}
        for activity in unique_activities_train:
            subset_train = df_train_emb[df_train_emb['Activity'] == activity]
            svgs_train = class_svgs_train[activity]
            if ID_column:
                id_values_train = list(self.df.loc[self.df.index.isin(subset_train.index)][ID_column])
                descriptions_train = [f'{ID_column}={value}' for value in id_values_train]  # No rounding for non-float IDs
            else:
                descriptions_train = [''] * len(subset_train)
            sources_train[activity] = ColumnDataSource(data=dict(
                x=subset_train['dim_1'],
                y=subset_train['dim_2'],
                desc=descriptions_train,
                svgs=svgs_train
            ))
    
        # Process the test set if it exists
        if test and self.df_test is not None:
            # Create SVGs for the test set
            svgs_test = [self.moltosvg(mol).data for mol in mols_test]
            
            # Create ColumnDataSource for the test set
            if ID_column:
                id_values_test = list(self.df_test[ID_column])
                descriptions_test = [f'{ID_column}={value}' for value in id_values_test]  # No rounding for non-float IDs
            else:
                descriptions_test = [''] * len(df_test_emb)
            source_test = ColumnDataSource(data=dict(
                x=df_test_emb['dim_1'],
                y=df_test_emb['dim_2'],
                desc=descriptions_test,
                svgs=svgs_test
            ))

        # Colors for training set class
        if len(unique_activities_train) == 2:
            colors_train = ['#00FF00', '#FF0000']  # Verde y Rojo
        else:
            colors_train = self.generate_colors(len(unique_activities_train))

        hover = HoverTool(tooltips="""
            <div>
                <div>@svgs{safe}</div>
                <div><span style="font-size: 17px; font-weight: bold;">@desc</span></div>
            </div>
        """)

        interactive_map = figure(
            width=800,
            height=800,
            tools=['reset', 'box_zoom', 'wheel_zoom', 'zoom_in', 'zoom_out', 'pan', hover],
            x_axis_label='Dimension 1',
            y_axis_label='Dimension 2'
        )

        # training set points
        for activity, color in zip(unique_activities_train, colors_train):
            interactive_map.circle('x', 'y', size=15, source=sources_train[activity], fill_alpha=0.2, color=color, legend_label=f'Class {activity}')

        # plot test set if exhist
        if test and self.df_test is not None:
            test_color = '#0000FF'  # Azul para el conjunto de test
            interactive_map.circle('x', 'y', size=15, source=source_test, fill_alpha=0.2, color=test_color, legend_label='Test set')

        interactive_map.legend.location = "top_right"
        interactive_map.legend.click_policy = "hide"

        if Save:
            output_file(path)
            save(interactive_map)
        else:
            show(interactive_map)