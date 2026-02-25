# ⚛️ Q-DISCOVER: Quantum-Inspired Ligand Discovery Platform

**Q-DISCOVER** is an enterprise-grade, full-stack chemoinformatics platform built to demonstrate the immense computational and economic advantages of Quantum Computing in pharmaceutical drug discovery. 

By mapping molecular databases to high-dimensional Hilbert spaces, Q-DISCOVER leverages **Qiskit** and **Grover's Amplitude Amplification** to reduce the algorithmic search complexity of finding viable drug candidates from Classical $O(N)$ to Quantum $O(\sqrt{N})$.

---

## 🚀 The Problem & The Quantum Solution
**The Classical Bottleneck:** Traditional computational screening evaluates molecules sequentially. For massive chemical databases, this linear $O(N)$ time complexity costs pharmaceutical companies billions of dollars in supercomputer time and years of R&D.

**The Quantum Advantage:** Q-DISCOVER utilizes IBM's Qiskit to create a uniform superposition of molecular states. By applying a multi-dimensional Oracle phase-flip and a Grover Diffuser, it forces destructive interference on invalid molecules and constructive interference on the target molecule, extracting the optimal drug candidate in a fraction of the time.

---

## 🧬 Key Platform Modules
1. **🏎️ Algorithm Race:** A live execution terminal visually pitting Classical sequential loops against Quantum matrix evolution.
2. **💰 Cost Simulator:** An interactive economic model proving that even if Quantum compute is 10x more expensive per step, the $O(\sqrt{N})$ scaling saves billions at scale.
3. **🌌 Hilbert Space:** 3D WebGL visualizations mapping the statevector probability distribution of the database.
4. **🛡️ Lipinski Filter & Oracle:** Enforces the "Rule of 5" to ensure extracted candidates are physically viable for human oral absorption.
5. **🧬 Tanimoto Similarity:** Utilizes RDKit Morgan Fingerprints (RAM-cached for 0.01s latency) to find structural alternative backups for top candidates.
6. **✍️ Molecule Designer:** A live Ketcher canvas allowing chemists to physically modify the extracted quantum hit and run real-time diagnostics against Oracle bounds.
7. **🎓 Student Tour Mode:** An interactive UI overlay that guides early researchers through the complex quantum mechanics and chemistry metrics.

---

## 🛠️ Technology Stack
* **Backend:** Python, FastAPI, Uvicorn
* **Quantum Engine:** IBM Qiskit (`QuantumCircuit`, `GroverOperator`, `Statevector`)
* **Chemoinformatics:** RDKit, Ketcher
* **Frontend:** HTML5, CSS3, Vanilla JavaScript, Plotly.js (Data Visualization), 3Dmol.js (3D Rendering)

> **Note on Hardware Optimization:** To ensure the live demonstration runs smoothly on classical laptop RAM without thermal-throttling, the active Qiskit simulation is scaled to a 10-Qubit subspace (1,000 molecules). 

---

## ⚙️ How to Run Locally

### 1. Install Dependencies
Ensure you have Python 3.8+ installed. Run the following command in your terminal:
```bash
pip install fastapi uvicorn pandas numpy rdkit qiskit ipython_genutils