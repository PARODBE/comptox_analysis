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
    
    def __init__(self,path:str=None):

        '''
        This library allows to draw the chemical space for a specific dataset, providing structure chemical information.

        ----

        Example of use:
        
        chemical_space_plotter('cox2_train.sdf').visualizer(type='morgan',ID_column='pchembl_value')
        
        '''
        
        self.df=PandasTools.LoadSDF(path,smilesName='Smiles',includeFingerprints=True,molColName='Structure')

    def _prepareMol(self, mol, kekulize):

        '''
        The objective of this function consists of preparing a molecule for visual representation by performing specific preprocessing steps
        '''
        
        mc = Chem.Mol(mol.ToBinary())
        if kekulize:
            try:
                Chem.Kekulize(mc)
            except:
                mc = Chem.Mol(mol.ToBinary())
        if not mc.GetNumConformers():
            rdDepictor.Compute2DCoords(mc)
        return mc

    def moltosvg(self, mol,molSize=(450,200),kekulize=True,drawer=None,**kwargs):

        '''
        The previous function is applied in moltosvg for visual representation of the different compounds.
        '''
        
        mc = self._prepareMol(mol,kekulize)
        if drawer is None:
            drawer = rdMolDraw2D.MolDraw2DSVG(molSize[0],molSize[1])
        drawer.DrawMolecule(mc,**kwargs)
        drawer.FinishDrawing()
        svg = drawer.GetDrawingText()
        
        return SVG(svg.replace('svg:',''))

    def generate_colors(self, num_colors):

        '''
        This function generates different colors for the different classes.
        '''
        
        palette = Category20[20]  # Usar una paleta con 20 colores
        color_cycle = cycle(palette)
        colors = [next(color_cycle) for _ in range(num_colors)]
        
        return colors

    def visualizer(self,type = 'morgan', hyp_fp : Dict = None, activity_column='Activity', ID_column=None):

        '''
        With this function you can visualize the chemical space of the selected dataframe. To complete this task you must select the type of
        descriptors --> molecular_descriptor (with the 208 rdkit variables), morgan and rdkit fingerprints. For the fingerprints you can specify the 
        different hyperparameters.
        Also, you must select the name of your activity_column and finally, an ID_column.
        '''

        mols=[Chem.MolFromSmiles(m) for m in self.df.Smiles]
        y=self.df[activity_column]

        if type=='molecular_descriptor':
            
            df_molecular_descriptors=pd.DataFrame([Descriptors.CalcMolDescriptors(m) for m in mols])
            sc=StandardScaler().fit(df_molecular_descriptors)
            X_desc_sc=sc.transform(df_molecular_descriptors)
            emb=PCA(n_components=2,random_state=46).fit_transform(X_desc_sc)

        elif type=='morgan':

            X_fp_morgan=np.array([AllChem.GetMorganFingerprintAsBitVect(mol,**hyp_fp) for mol in 
                                  [Chem.MolFromSmiles(m) for m in list(self.df.Smiles)]])
            emb = PredictableTSNE(estimator=MLPRegressor(random_state=46), transformer=TSNE(random_state=46) ).fit_transform(X_fp_morgan,y)
            
        elif type=='rdkit':

            X_fp_rdkit = np.array([AllChem.RDKFingerprint(mol,**hyp_fp) for mol in [Chem.MolFromSmiles(m) for m in list(self.df.Smiles)]])
            emb = PredictableTSNE(estimator=MLPRegressor(random_state=46), transformer=TSNE(random_state=46) ).fit_transform(X_fp_rdkit,y)

            
        df_emb = pd.DataFrame(np.c_[emb, y], index=self.df.index.values, columns=['dim_1', 'dim_2', 'Activity'])
        unique_activities = y.unique()
        class_mols = {activity: [Chem.MolFromSmiles(smiles) for smiles in self.df[self.df[activity_column] == activity]['Smiles']] for activity in unique_activities}
        class_svgs = {activity: [self.moltosvg(mol).data for mol in mols] for activity, mols in class_mols.items()}
        
        sources = {}
        for activity in unique_activities:
            subset = df_emb[df_emb['Activity'] == activity]
            svgs = class_svgs[activity]
            id_values = list(self.df.loc[self.df.index.isin(subset.index)][ID_column])
            descriptions = [f'{ID_column}={value}' for value in id_values]
            sources[activity] = ColumnDataSource(data=dict(
                x=subset['dim_1'],
                y=subset['dim_2'],
                desc=descriptions,
                svgs=svgs
            ))

        if len(unique_activities) == 2:
            colors = ['#00FF00', '#FF0000']  # Verde y Rojo
        else:
            colors = self.generate_colors(len(unique_activities))

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

        for activity, color in zip(unique_activities, colors):
            interactive_map.circle('x', 'y', size=15, source=sources[activity], fill_alpha=0.2, color=color, legend_label=f'Class {activity}')

        interactive_map.legend.location = "top_right"
        interactive_map.legend.click_policy = "hide"

        show(interactive_map)