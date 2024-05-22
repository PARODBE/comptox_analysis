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
    
    def __init__(self, smile:List[str]=None):
        self.smile = smile
        self.mol = [Chem.MolFromSmiles(sm) for sm in self.smile]
        self.bi_list = []
        self.info_ex_list = []
        self.smarts_list = []
        self.hit_ats_list_list = []
        self.hit_bonds_list_list = []
    
    def render_image(self, size:Tuple[int,int]=(600,600), indexes:bool=True, number:int=None) -> Image.Image:
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
        d2d.DrawMolecule(self.mol[number])
        d2d.FinishDrawing()
        return SVG(d2d.GetDrawingText())

    def fragmentation(self, n: int = None, number: int = None) -> Tuple[Image.Image, Dict]:
        '''
        This function gives your molecular fragmentation
        ----
        Hyperparameters:
    
        n = Number of fragments
        '''
        # if len(self.mol)==1:

        tpls_list = []
        for mol in self.mol:
            bi = {}
            fp = rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, radius=2, bitInfo=bi, useFeatures=True)
            # show 10 of the set bits:
            tpls = [(mol, x, bi) for x in fp.GetOnBits()]
            self.bi_list.append(bi)
            tpls_list.append(tpls)
            
        return (Draw.DrawMorganBits(tpls_list[number][:n], molsPerRow=4, legends=[str(x) for x in fp.GetOnBits()][:n]), self.bi_list[number])

    def getSubstructSmi(self, mol, atomID, radius):
        '''
        This function gives your text information about molecular fragmentation
        
        Source code : https://rdkit.blogspot.com/2016/03/explaining-morgan-similarity.html
        '''

        if radius > 0:
            env = Chem.FindAtomEnvironmentOfRadiusN(mol, radius, atomID)
            atomsToUse = []
            for b in env:
                atomsToUse.append(mol.GetBondWithIdx(b).GetBeginAtomIdx())
                atomsToUse.append(mol.GetBondWithIdx(b).GetEndAtomIdx())
            atomsToUse = list(set(atomsToUse))
        else:
            atomsToUse = [atomID]
            env = None
        symbols = []
        for atom in mol.GetAtoms():
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
        
        smi = Chem.MolFragmentToSmiles(mol, atomsToUse, bondsToUse=env, allHsExplicit=True, allBondsExplicit=True, rootedAtAtom=atomID)
        smi2 = Chem.MolFragmentToSmiles(mol, atomsToUse, bondsToUse=env, atomSymbols=symbols, allBondsExplicit=True, rootedAtAtom=atomID)
        
        return smi, smi2
    
    def highlighting(self, fingerprint_numbers: List[List[int]] = None) -> Image.Image:
        '''
        This function draws the selected fingerprint(s) for the selected compound(s)
        '''
        for mol in self.mol:
            bi = {}
            fp = rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, radius=2, bitInfo=bi, useFeatures=True)
            # show 10 of the set bits:
            self.bi_list.append(bi)
            
        for bi, mol in zip(self.bi_list, self.mol):
            info_ex = []
            for bitId, atoms in bi.items():
                exampleAtom, exampleRadius = atoms[0]
                description = self.getSubstructSmi(mol, exampleAtom, exampleRadius)
                info_ex.append((bitId, exampleRadius, description[0], description[1]))
            self.info_ex_list.append(info_ex)

        for info_ex, fingerprint_number in zip(self.info_ex_list, fingerprint_numbers):
            smarts = []
            for finger in fingerprint_number:
                for i in info_ex:
                    if i[0] == finger:
                        smart = i[2]
                        smarts.append(smart)
            self.smarts_list.append(smarts)

        for mol, smarts in zip(self.mol, self.smarts_list):
            hit_ats_list = []
            hit_bonds_list = []
            for smart in smarts:
                patt = Chem.MolFromSmarts(smart)
                hit_ats = list(mol.GetSubstructMatch(patt))
                hit_bonds = []
                for bond in patt.GetBonds():
                    aid1 = hit_ats[bond.GetBeginAtomIdx()]
                    aid2 = hit_ats[bond.GetEndAtomIdx()]
                    hit_bonds.append(mol.GetBondBetweenAtoms(aid1, aid2).GetIdx())
                hit_ats_list.append(hit_ats)
                hit_bonds_list.append(hit_bonds)
            self.hit_ats_list_list.append(hit_ats_list)
            self.hit_bonds_list_list.append(hit_bonds_list)

        if len(self.mol) == 1:
            result_dict_ats = defaultdict(list)
            result_dict_bonds = defaultdict(list)
            
            default_color = (1.0, 0.0, 0.0, 0.4)
            for hit_ats, hit_bonds in zip(self.hit_ats_list_list[0], self.hit_bonds_list_list[0]):
                for bit in hit_ats:
                    if bit not in result_dict_ats:
                        result_dict_ats[bit].append(default_color)
                for bit in hit_bonds:
                    if bit not in result_dict_bonds:
                        result_dict_bonds[bit].append(default_color)

            result_ats = dict(result_dict_ats)
            result_bonds = dict(result_dict_bonds)
            
            d2d = rdMolDraw2D.MolDraw2DCairo(600, 600)
            d2d.drawOptions().addAtomIndices = True
            d2d.DrawMoleculeWithHighlights(self.mol[0], "", result_ats, result_bonds, {}, {})
            d2d.FinishDrawing()
            bio = io.BytesIO(d2d.GetDrawingText())
            return Image.open(bio)
        else:
            default_color = (1.0, 0.0, 0.0, 0.4)
            colors = [default_color] * len(self.mol)
            
            highlightAtomLists = []
            highlightBondLists = []
            for hit_ats_list, hit_bonds_list, color in zip(self.hit_ats_list_list, self.hit_bonds_list_list, colors):
                result_dict_ats = defaultdict(list)
                result_dict_bonds = defaultdict(list)
                for hit_ats, hit_bonds in zip(hit_ats_list, hit_bonds_list):
                    for bit in hit_ats:
                        result_dict_ats[bit].append(color)
                    for bit in hit_bonds:
                        result_dict_bonds[bit].append(color)
                highlightAtomLists.append(dict(result_dict_ats))
                highlightBondLists.append(dict(result_dict_bonds))
            
            return Draw.MolsToGridImage(self.mol, highlightAtomLists=highlightAtomLists, highlightBondLists=highlightBondLists, subImgSize=(600,600))
