// Relative path to seamlessly connect to FastAPI on the exact same port
const API_URL = "/api";

document.addEventListener("DOMContentLoaded", () => {
    fetchMolecules();
    document.getElementById("run-sim-btn").addEventListener("click", runQuantumSimulation);
});

async function fetchMolecules() {
    try {
        const response = await fetch(`${API_URL}/molecules`);
        const data = await response.json();
        updateTable(data);
        renderPlot(data);
    } catch (error) {
        console.error("Error fetching data:", error);
    }
}

async function runQuantumSimulation() {
    const maxMw = document.getElementById("max-mw").value;
    const groverSteps = document.getElementById("grover-steps").value;
    const statusMsg = document.getElementById("status-msg");

    statusMsg.innerText = "Initializing Qiskit Statevector Engine...";

    try {
        const response = await fetch(`${API_URL}/simulate`, {
            method: "POST",
            headers: { "Content-Type": "application/json" },
            body: JSON.stringify({ max_mw: maxMw, grover_steps: groverSteps })
        });

        const result = await response.json();
        statusMsg.innerText = result.message;

        updateTable(result.data);
        renderPlot(result.data);
    } catch (error) {
        statusMsg.innerText = "Simulation failed. Ensure Uvicorn is running.";
        console.error(error);
    }
}

function updateTable(data) {
    const tbody = document.getElementById("table-body");
    tbody.innerHTML = "";

    data.sort((a, b) => b.probability - a.probability);

    data.forEach(mol => {
        const tr = document.createElement("tr");
        const statusBadge = mol.status === "PASS"
            ? `<span class="badge-pass">COMPLIANT</span>`
            : `<span class="badge-fail">FAIL (${mol.violations})</span>`;

        tr.innerHTML = `
            <td><strong>${mol.name}</strong></td>
            <td>${mol.mw}</td>
            <td>${mol.logp}</td>
            <td>${mol.hbd}</td>
            <td>${mol.hba}</td>
            <td>${statusBadge}</td>
            <td style="font-family: monospace; color: #f59e0b;">
                ${(mol.probability * 100).toFixed(4)}%
            </td>
        `;
        tbody.appendChild(tr);
    });
}

function renderPlot(data) {
    const x = data.map(d => d.mw);
    const y = data.map(d => d.logp);
    const z = data.map(d => d.probability);
    const texts = data.map(d => d.name);

    const sizes = data.map(d => Math.max(d.probability * 100, 6));

    const colors = data.map(d => {
        if (d.classification === '★ MEASURED WAVE COLLAPSE') return '#f59e0b';
        if (d.classification === 'Amplified Target') return '#3b82f6';
        return '#475569';
    });

    const trace = {
        x: x, y: y, z: z,
        mode: 'markers',
        marker: { size: sizes, color: colors, opacity: 0.9, line: { width: 0 } },
        text: texts, type: 'scatter3d', hoverinfo: 'text'
    };

    const layout = {
        margin: { l: 0, r: 0, b: 0, t: 0 },
        paper_bgcolor: 'rgba(0,0,0,0)',
        scene: {
            xaxis: { title: 'Mass (Da)', backgroundcolor: 'rgba(100,100,100,0.05)', color: 'white' },
            yaxis: { title: 'LogP', backgroundcolor: 'rgba(100,100,100,0.05)', color: 'white' },
            zaxis: { title: '|Ψ|² Prob', backgroundcolor: 'rgba(100,100,100,0.05)', color: 'white' },
            bgcolor: 'rgba(0,0,0,0)'
        }
    };

    Plotly.newPlot('plot3d', [trace], layout, {displayModeBar: false});
}