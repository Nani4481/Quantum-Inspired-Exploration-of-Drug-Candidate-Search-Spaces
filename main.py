import streamlit as st
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors

# ==========================================
# 0. CONFIG & THEME
# ==========================================
st.set_page_config(page_title="Q-DISCOVER | V1 MVP", layout="wide", page_icon="⚛️")
st.title("⚛️ Q-DISCOVER - Initial Prototype")
st.write("Version 1: Basic Chemoinformatics Validation")

# ==========================================
# 1. BASE DATASET (10 MOLECULES)
# ==========================================
sample_drugs = {
    "Aspirin": "CC(=O)OC1=CC=CC=C1C(=O)O",
    "Paracetamol": "CC(=O)NC1=CC=C(O)C=C1",
    "Ibuprofen": "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",
    "Imatinib": "CC1=C(C=C(C=C1)NC(=O)C2=CC=C(C=C2)CN3CCN(CC3)C)NC4=NC=CC(=N4)C5=CN=CC=C5",
    "Fluoxetine": "CNCCC(C1=CC=CC=C1)OC2=CC=C(C=C2)C(F)(F)F",
    "Omeprazole": "CC1=CN=C(C(=C1OC)C)CS(=O)C2=NC3=C(N2)C=C(C=C3)OC",
    "Amoxicillin": "CC1(C(N2C(S1)C(C2=O)NC(=O)C(C3=CC=C(C=C3)O)N)C(=O)O)C",
    "Sildenafil": "CCCC1=NN(C2=C1N=C(NC2=O)C3=C(C=CC(=C3)S(=O)(=O)N4CCN(CC4)C)OCC)C",
    "Metformin": "CN(C)C(=N)NC(=N)N",
    "Caffeine": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
}

# Process SMILES strings through RDKit
data = []
for name, smi in sample_drugs.items():
    mol = Chem.MolFromSmiles(smi)
    if mol:
        data.append({
            "Name": name,
            "SMILES": smi,
            "Molecular Weight (Da)": Descriptors.MolWt(mol),
            "LogP (Lipophilicity)": Descriptors.MolLogP(mol)
        })

df = pd.DataFrame(data)

# ==========================================
# 2. UI RENDERING
# ==========================================
st.markdown("### Processed Molecular Search Space")
st.dataframe(df, use_container_width=True)