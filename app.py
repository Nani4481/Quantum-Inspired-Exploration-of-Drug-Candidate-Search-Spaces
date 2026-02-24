import os
import uvicorn
from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from fastapi.staticfiles import StaticFiles
from fastapi.responses import FileResponse
from pydantic import BaseModel
import pandas as pd
import numpy as np
import math

from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

from qiskit import QuantumCircuit
from qiskit.circuit.library import GroverOperator, DiagonalGate
from qiskit.quantum_info import Statevector

app = FastAPI(title="Q-DISCOVER API")

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

sample_drugs = {
    "Aspirin": "CC(=O)OC1=CC=CC=C1C(=O)O", "Paracetamol": "CC(=O)NC1=CC=C(O)C=C1",
    "Ibuprofen": "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",
    "Imatinib": "CC1=C(C=C(C=C1)NC(=O)C2=CC=C(C=C2)CN3CCN(CC3)C)NC4=NC=CC(=N4)C5=CN=CC=C5",
    "Fluoxetine": "CNCCC(C1=CC=CC=C1)OC2=CC=C(C=C2)C(F)(F)F",
    "Omeprazole": "CC1=CN=C(C(=C1OC)C)CS(=O)C2=NC3=C(N2)C=C(C=C3)OC",
    "Amoxicillin": "CC1(C(N2C(S1)C(C2=O)NC(=O)C(C3=CC=C(C=C3)O)N)C(=O)O)C",
    "Sildenafil": "CCCC1=NN(C2=C1N=C(NC2=O)C3=C(C=CC(=C3)S(=O)(=O)N4CCN(CC4)C)OCC)C",
    "Metformin": "CN(C)C(=N)NC(=N)N", "Caffeine": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
}


def get_base_data():
    data = []
    for name, smi in sample_drugs.items():
        mol = Chem.MolFromSmiles(smi)
        if mol:
            mw = Descriptors.MolWt(mol)
            logp = Descriptors.MolLogP(mol)
            hbd = rdMolDescriptors.CalcNumHBD(mol)
            hba = rdMolDescriptors.CalcNumHBA(mol)

            violations = sum([mw > 500, logp > 5, hbd > 5, hba > 10])
            lipinski_status = "PASS" if violations <= 1 else "FAIL"

            data.append({
                "name": name, "smiles": smi, "mw": round(mw, 2), "logp": round(logp, 2),
                "hbd": hbd, "hba": hba, "violations": violations, "status": lipinski_status,
                "probability": 0.0
            })
    return pd.DataFrame(data)


class SimulationRequest(BaseModel):
    max_mw: float = 300.0
    grover_steps: int = 2


@app.get("/api/molecules")
def get_molecules():
    df = get_base_data()
    return df.to_dict(orient='records')


@app.post("/api/simulate")
def run_quantum_simulation(req: SimulationRequest):
    df = get_base_data()
    N = len(df)

    oracle_mask = (df['mw'] <= req.max_mw).values
    targets_count = np.sum(oracle_mask)

    if targets_count == 0:
        return {"message": "No targets matched.", "data": df.to_dict(orient='records')}

    num_qubits = math.ceil(math.log2(N))
    dim = 2 ** num_qubits
    phases = np.ones(dim, dtype=np.complex128)

    targets = np.where(oracle_mask)[0]
    for t in targets:
        phases[t] = -1.0

    oracle_gate = DiagonalGate(phases.tolist())
    oracle_qc = QuantumCircuit(num_qubits)
    oracle_qc.append(oracle_gate, range(num_qubits))
    grover_op = GroverOperator(oracle_qc)

    qc = QuantumCircuit(num_qubits)
    qc.h(range(num_qubits))
    sv = Statevector(qc)

    for _ in range(req.grover_steps):
        sv = sv.evolve(grover_op)

    df['probability'] = sv.probabilities()[:N]

    df['classification'] = 'Suppressed Noise'
    uniform_prob = 1.0 / dim
    df.loc[df['probability'] > uniform_prob * 1.5, 'classification'] = 'Amplified Target'
    top_hit_idx = df['probability'].idxmax()
    df.loc[top_hit_idx, 'classification'] = '★ MEASURED WAVE COLLAPSE'

    return {
        "message": f"Wavefunction Collapsed! Found {targets_count} targets.",
        "data": df.to_dict(orient='records')
    }


# ==========================================
# EXPLICIT FRONTEND ROUTING FIX
# ==========================================
# ... (all your API routes remain the same) ...

# Final Corrected Path Logic
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
# Points to the folder containing your .css and .js files
STATIC_FILES_DIR = os.path.join(BASE_DIR, "frontend", "static")

# Mount the folder
app.mount("/static", StaticFiles(directory=STATIC_FILES_DIR), name="static")

@app.get("/")
def serve_index():
    # Points to index.html which is inside the 'frontend' folder
    return FileResponse(os.path.join(BASE_DIR, "frontend", "index.html"))

if __name__ == '__main__':
    uvicorn.run("app:app", host="127.0.0.1", port=8000, reload=True)