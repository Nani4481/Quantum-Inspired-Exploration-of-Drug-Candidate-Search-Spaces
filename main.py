from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel
import pandas as pd
import numpy as np
import math
import os
import warnings
import uvicorn

warnings.filterwarnings("ignore", category=DeprecationWarning)

from rdkit import Chem
from rdkit.Chem import Descriptors, QED, AllChem, rdMolDescriptors
from rdkit.Chem import DataStructs
from qiskit import QuantumCircuit
from qiskit.circuit.library import GroverOperator, DiagonalGate
from qiskit.quantum_info import Statevector

app = FastAPI(title="Q-DISCOVER Enterprise Backend")

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_methods=["*"],
    allow_headers=["*"],
)

GLOBAL_DF = pd.DataFrame()

@app.on_event("startup")
def load_data():
    global GLOBAL_DF
    data_file = "real_10k_dataset.csv"

    real_drugs = {
        "Aspirin": "CC(=O)OC1=CC=CC=C1C(=O)O", "Paracetamol": "CC(=O)NC1=CC=C(O)C=C1",
        "Ibuprofen": "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",
        "Imatinib": "CC1=C(C=C(C=C1)NC(=O)C2=CC=C(C=C2)CN3CCN(CC3)C)NC4=NC=CC(=N4)C5=CN=CC=C5",
        "Fluoxetine": "CNCCC(C1=CC=CC=C1)OC2=CC=C(C=C2)C(F)(F)F",
        "Omeprazole": "CC1=CN=C(C(=C1OC)C)CS(=O)C2=NC3=C(N2)C=C(C=C3)OC",
        "Amoxicillin": "CC1(C(N2C(S1)C(C2=O)NC(=O)C(C3=CC=C(C=C3)O)N)C(=O)O)C",
        "Sildenafil": "CCCC1=NN(C2=C1N=C(NC2=O)C3=C(C=CC(=C3)S(=O)(=O)N4CCN(CC4)C)OCC)C",
        "Metformin": "CN(C)C(=N)NC(=N)N",
        "Atorvastatin": "CC(C)C1=C(C(=C(N1CCC(CC(CC(=O)O)O)O)C2=CC=C(C=C2)F)C3=CC=CC=C3)C(=O)NC4=CC=CC=C4",
        "Caffeine": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C", "Lisinopril": "OC(=O)C(CCCCN)NC(C(=O)N1CCCC1C(=O)O)CC2=CC=CC=C2",
    }

    base_data = []
    for name, smi in real_drugs.items():
        mol = Chem.MolFromSmiles(smi)
        if mol:
            base_data.append({
                'Name': name, 'SMILES': smi, 'Is_Real': True, 'MW': Descriptors.MolWt(mol),
                'LogP': Descriptors.MolLogP(mol), 'QED': QED.qed(mol), 'TPSA': Descriptors.TPSA(mol),
                'HBD': rdMolDescriptors.CalcNumHBD(mol), 'HBA': rdMolDescriptors.CalcNumHBA(mol),
                'RotBonds': rdMolDescriptors.CalcNumRotatableBonds(mol), 'Rings': rdMolDescriptors.CalcNumRings(mol),
            })

    if os.path.exists(data_file):
        df_real = pd.read_csv(data_file)
        GLOBAL_DF = pd.concat([pd.DataFrame(base_data), df_real], ignore_index=True).drop_duplicates(subset=['SMILES']).reset_index(drop=True)
        print(f"✅ Loaded {len(GLOBAL_DF)} molecules from local cache.")
        return

    print("⏳ Downloading dataset from AWS. This may take a minute...")
    try:
        df_raw = pd.read_csv("https://deepchemdata.s3-us-west-1.amazonaws.com/datasets/HIV.csv")
        smiles_list = df_raw['smiles'].dropna().unique()[:10000]
        real_data = []
        for i, smi in enumerate(smiles_list):
            mol = Chem.MolFromSmiles(smi)
            if mol:
                try:
                    real_data.append({
                        'Name': f"ChemBase-{i}", 'SMILES': smi, 'Is_Real': True,
                        'MW': Descriptors.MolWt(mol), 'LogP': Descriptors.MolLogP(mol),
                        'QED': QED.qed(mol), 'TPSA': Descriptors.TPSA(mol),
                        'HBD': rdMolDescriptors.CalcNumHBD(mol), 'HBA': rdMolDescriptors.CalcNumHBA(mol),
                        'RotBonds': rdMolDescriptors.CalcNumRotatableBonds(mol),
                        'Rings': rdMolDescriptors.CalcNumRings(mol),
                    })
                except: continue
        df_fetched = pd.DataFrame(real_data)
        df_fetched.to_csv(data_file, index=False)
        GLOBAL_DF = pd.concat([pd.DataFrame(base_data), df_fetched], ignore_index=True).drop_duplicates(subset=['SMILES']).reset_index(drop=True)
        print(f"✅ Successfully downloaded {len(GLOBAL_DF)} molecules.")
    except Exception as e:
        print(f"⚠️ AWS Download failed ({e}). Generating synthetic vector space...")
        real_data = []
        np.random.seed(42)
        for i in range(9988):
            real_data.append({
                'Name': f"SYN-{(i + 1):05d}", 'SMILES': "C", 'Is_Real': False,
                'MW': abs(np.random.normal(550, 200)), 'LogP': np.random.normal(5.5, 3.0),
                'QED': min(max(np.random.beta(1.5, 8.0), 0.01), 0.99), 'TPSA': abs(np.random.normal(120, 60)),
                'HBD': int(np.random.poisson(2)), 'HBA': int(np.random.poisson(4)),
                'RotBonds': int(np.random.poisson(5)), 'Rings': int(np.random.poisson(2)),
            })
        GLOBAL_DF = pd.concat([pd.DataFrame(base_data), pd.DataFrame(real_data)], ignore_index=True).drop_duplicates(subset=['SMILES']).reset_index(drop=True)

class OracleConfig(BaseModel):
    t_qed: float
    min_mw: float
    max_mw: float
    min_logp: float
    max_logp: float
    min_tpsa: float
    max_tpsa: float
    max_hbd: int
    max_hba: int
    max_rotbonds: int
    max_rings: int
    grover_steps: int
    noise_rate: float

@app.post("/api/simulate")
def run_quantum_simulation(config: OracleConfig):
    global GLOBAL_DF
    local_df = GLOBAL_DF.copy()
    N = len(local_df)

    range_mask = (
            (local_df['QED'] >= config.t_qed) & (local_df['MW'] >= config.min_mw) & (local_df['MW'] <= config.max_mw) &
            (local_df['LogP'] >= config.min_logp) & (local_df['LogP'] <= config.max_logp) &
            (local_df['TPSA'] >= config.min_tpsa) & (local_df['TPSA'] <= config.max_tpsa) &
            (local_df['HBD'] <= config.max_hbd) & (local_df['HBA'] <= config.max_hba) &
            (local_df['RotBonds'] <= config.max_rotbonds) & (local_df['Rings'] <= config.max_rings)
    )
    oracle_mask = range_mask.values
    targets_count = int(np.sum(oracle_mask))

    num_qubits = math.ceil(math.log2(N))
    dim = 2 ** num_qubits
    phases = np.ones(dim, dtype=np.complex128)

    if targets_count > 0:
        targets = np.where(oracle_mask)[0]
        valid_qeds = local_df.loc[targets, 'QED'].values
        if len(valid_qeds) > 0:
            max_qed = valid_qeds.max()
            norm_scores = valid_qeds / max_qed
            for t, score in zip(targets, norm_scores):
                phases[t] = np.exp(1j * np.pi * score)

    oracle_gate = DiagonalGate(phases.tolist())
    oracle_qc = QuantumCircuit(num_qubits)
    oracle_qc.append(oracle_gate, range(num_qubits))

    grover_op = GroverOperator(oracle_qc)
    qc = QuantumCircuit(num_qubits)
    qc.h(range(num_qubits))
    sv = Statevector(qc)

    psi_history = [sv.probabilities()[:N].tolist()]

    for step in range(config.grover_steps):
        sv = sv.evolve(grover_op)
        probs = sv.probabilities()[:N]
        if config.noise_rate > 0:
            probs = (probs * (1.0 - config.noise_rate)) + (np.ones(N) * (1.0 / N) * config.noise_rate)
        psi_history.append(probs.tolist())

    local_df['Probability'] = psi_history[-1]
    local_df = local_df.fillna(0).replace([np.inf, -np.inf], 0)

    top_hit = local_df.iloc[int(local_df['Probability'].idxmax())]
    M_approx = max(1, targets_count)
    opt_steps = max(1, int((np.pi / 4) * math.sqrt(dim / M_approx)))

    uniform_prob = 1.0 / N
    local_df['Classification'] = 'Suppressed Noise'
    local_df.loc[local_df['Probability'] > uniform_prob * 1.5, 'Classification'] = 'Amplified Target(s)'
    local_df.loc[local_df['Name'] == top_hit['Name'], 'Classification'] = '★ MEASURED WAVE COLLAPSE'

    return {
        "N": N, "dim": dim, "num_qubits": num_qubits, "targets_count": targets_count,
        "opt_steps": opt_steps, "top_hit": top_hit.to_dict(), "psi_history": psi_history,
        "display_df": local_df.sort_values('Probability', ascending=False).head(1000).to_dict(orient="records")
    }

class SmilesQuery(BaseModel): smiles: str

@app.post("/api/similarity")
def similarity_search(query: SmilesQuery):
    query_mol = Chem.MolFromSmiles(query.smiles)
    if not query_mol: raise HTTPException(status_code=400, detail="Invalid SMILES")
    query_fp = AllChem.GetMorganFingerprintAsBitVect(query_mol, 2, nBits=1024)

    real_drugs = GLOBAL_DF[GLOBAL_DF['Is_Real']].copy()
    sims = []
    for _, row in real_drugs.iterrows():
        mol_b = Chem.MolFromSmiles(row['SMILES'])
        if mol_b:
            score = DataStructs.TanimotoSimilarity(query_fp, AllChem.GetMorganFingerprintAsBitVect(mol_b, 2, nBits=1024))
            sims.append({'Name': row['Name'], 'SMILES': row['SMILES'], 'Tanimoto': round(score, 4)})

    matches = pd.DataFrame(sims).sort_values('Tanimoto', ascending=False).head(10).fillna(0)
    return {"matches": matches.to_dict(orient="records")}

@app.get("/api/3dmol")
def get_3dmol(smiles: str):
    try:
        mol = Chem.MolFromSmiles(smiles)
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol)
        AllChem.MMFFOptimizeMolecule(mol)
        return {"mol_block": Chem.MolToMolBlock(mol)}
    except:
        raise HTTPException(status_code=400, detail="3D generation failed")

@app.post("/api/molecule-props")
def get_props(query: SmilesQuery):
    mol = Chem.MolFromSmiles(query.smiles)
    if not mol: raise HTTPException(status_code=400, detail="Invalid SMILES")
    Chem.SanitizeMol(mol)
    props = {
        'MW': float(Descriptors.MolWt(mol)), 'LogP': float(Descriptors.MolLogP(mol)),
        'QED': float(QED.qed(mol)), 'TPSA': float(Descriptors.TPSA(mol)),
        'HBD': int(rdMolDescriptors.CalcNumHBD(mol)), 'HBA': int(rdMolDescriptors.CalcNumHBA(mol)),
        'RotBonds': int(rdMolDescriptors.CalcNumRotatableBonds(mol)), 'Rings': int(rdMolDescriptors.CalcNumRings(mol)),
    }
    checks = {"MW <= 500 Da": props['MW'] <= 500, "LogP <= 5": props['LogP'] <= 5, "HBD <= 5": props['HBD'] <= 5, "HBA <= 10": props['HBA'] <= 10}
    return {"props": props, "checks": checks, "violations": sum(1 for v in checks.values() if not v)}

if __name__ == "__main__":
    print("\n⚛️  Starting Q-DISCOVER Fast API Backend...")
    uvicorn.run(app, host="127.0.0.1", port=8000)