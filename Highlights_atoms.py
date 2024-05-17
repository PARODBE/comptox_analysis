from rdkit.Chem import PandasTools, Descriptors
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import DataStructs
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import Draw
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from IPython.display import SVG
from PIL import Image
import io
from typing import Tuple, List, Dict
from collections import defaultdict
import warnings

warnings.filterwarnings("ignore")



class highlighter:
    
    def __init__(self, smile:str=None):
        
        self.smile=smile
        self.mol = Chem.MolFromSmiles(self.smile)
        self.bi = {}
        self.info_ex = []
        self.smarts = []
        self.hit_ats=[]
        self.hit_bonds = []
        self.hit_ats_list = []
        self.hit_bonds_list = []
    
    def render_image(self, size:Tuple[int,int]=(600,600), indexes:bool=True) -> Image.Image:
        '''
        This function render your image.
        You have to introduce the smile of the compound that you want to render.
        ----
        Hyperparameters:
        
        size = Structure size
        indexes = atom number
        '''
        
        width, height = size
        d2d = rdMolDraw2D.MolDraw2DSVG(width, height)
        d2d.drawOptions().addAtomIndices = indexes
        d2d.DrawMolecule(self.mol)
        d2d.FinishDrawing()
        return SVG(d2d.GetDrawingText())

    def fragmentation(self, n: int = None) -> Tuple[Image.Image, Dict]:
        '''
        This function gives your molecular fragmentation
        
        Hyperparameters:
    
        n = Number of fragments
        '''
        fp = rdMolDescriptors.GetMorganFingerprintAsBitVect(self.mol, radius=2, bitInfo=self.bi, useFeatures=True)
        # show 10 of the set bits:
        tpls = [(self.mol, x, self.bi) for x in fp.GetOnBits()]
        return (Draw.DrawMorganBits(tpls[:n], molsPerRow=4, legends=[str(x) for x in fp.GetOnBits()][:n]), self.bi)

    def getSubstructSmi(self,mol,atomID,radius):
        '''
        This function gives your text information about molecular fragmentation
        
        Source code : https://rdkit.blogspot.com/2016/03/explaining-morgan-similarity.html
        '''

        if radius > 0:
            env = Chem.FindAtomEnvironmentOfRadiusN(self.mol, radius, atomID)
            atomsToUse = []
            for b in env:
                atomsToUse.append(self.mol.GetBondWithIdx(b).GetBeginAtomIdx())
                atomsToUse.append(self.mol.GetBondWithIdx(b).GetEndAtomIdx())
            atomsToUse = list(set(atomsToUse))
        else:
            atomsToUse = [atomID]
            env = None
        symbols = []
        for atom in self.mol.GetAtoms():
            deg = atom.GetDegree()
            isInRing = atom.IsInRing()
            nHs = atom.GetTotalNumHs()
            symbol = '[' + atom.GetSmarts()
            if nHs: 
                symbol += 'H'
                if nHs > 1:
                    symbol += '%d' % nHs
            if isInRing:
                symbol += ';R'
            else:
                symbol += ';!R'
            symbol += ';D%d' % deg
            symbol += "]"
            symbols.append(symbol)
        
        smi = Chem.MolFragmentToSmiles(self.mol, atomsToUse, bondsToUse=env, allHsExplicit=True, allBondsExplicit=True, rootedAtAtom=atomID)
        smi2 = Chem.MolFragmentToSmiles(self.mol, atomsToUse, bondsToUse=env, atomSymbols=symbols, allBondsExplicit=True, rootedAtAtom=atomID)
        
        return smi, smi2

    def get_smart(self):

        '''
        This function allows to obtain fragmentation option accordint to a given fingerprint
        '''
        
        for bitId, atoms in self.bi.items():
            exampleAtom, exampleRadius = atoms[0]
            description = self.getSubstructSmi(self.mol, exampleAtom, exampleRadius)
            self.info_ex.append((bitId, exampleRadius, description[0], description[1]))

    def finger(self, fingerprint_numbers=None):

        '''
        The result of previous function
        '''
        
        self.smarts.clear()  # Borra todos los elementos existentes en la lista self.smarts
        for fingerprint_number in fingerprint_numbers:
            for i in self.info_ex:
                if i[0] == fingerprint_number:
                    smart = i[2]
                    self.smarts.append(smart)
        # return self.smarts

    def atom_bond(self):
        
        for smart in self.smarts:
            patt = Chem.MolFromSmarts(smart)
            self.hit_ats = list(self.mol.GetSubstructMatch(patt))
            for bond in patt.GetBonds():
                aid1 = self.hit_ats[bond.GetBeginAtomIdx()]
                aid2 = self.hit_ats[bond.GetEndAtomIdx()]
                self.hit_bonds.append(self.mol.GetBondBetweenAtoms(aid1, aid2).GetIdx())
            self.hit_ats_list.append(self.hit_ats)
            self.hit_bonds_list.append(self.hit_bonds)
    
        # return self.hit_ats_list, self.hit_bonds_list

    def highlighting(self, colors=None):
        result_dict_ats = defaultdict(list)
        result_dict_bonds = defaultdict(list)
    
        if colors:
            if len(self.hit_ats_list) == 1:
                color = colors[0]
                for bit in hit_ats_list[0]:
                    result_dict_ats[bit].append(color)
                for bit in hit_bonds_list[0]:
                    result_dict_bonds[bit].append(color)
            else:
                for i, (hit_ats, hit_bonds) in enumerate(zip(self.hit_ats_list, self.hit_bonds_list)):
                    color = colors[i % len(colors)]  # Usar el color correspondiente según el índice
                    for bit in hit_ats:
                        result_dict_ats[bit].append(color)
                    for bit in hit_bonds:
                        result_dict_bonds[bit].append(color)
        else:
            # Si no se proporcionan colores, usar el color por defecto (rojo)
            default_color = (1.0, 0.0, 0.0, 0.4)
            for hit_ats, hit_bonds in zip(self.hit_ats_list, self.hit_bonds_list):
                for bit in hit_ats:
                    if bit not in result_dict_ats:
                        result_dict_ats[bit].append(default_color)
                for bit in hit_bonds:
                    if bit not in result_dict_bonds:
                        result_dict_bonds[bit].append(default_color)
    
        # Convertir el defaultdict a dict
        result_ats = dict(result_dict_ats)
        result_bonds = dict(result_dict_bonds)
    
        # Crear la imagen resaltada
        d2d = rdMolDraw2D.MolDraw2DCairo(600, 600)
        d2d.drawOptions().addAtomIndices = True
        d2d.DrawMoleculeWithHighlights(self.mol, "", result_ats, result_bonds, {}, {})
        d2d.FinishDrawing()
        bio = io.BytesIO(d2d.GetDrawingText())
        return Image.open(bio)