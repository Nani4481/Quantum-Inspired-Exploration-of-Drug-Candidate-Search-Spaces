const API_URL = "http://127.0.0.1:8000/api";
let globalSimData = null;
let currentActivePage = "Home";
let costAnimationInterval = null;
let lastLoadedSmiles = null;
let studentModeActive = false;

// ==========================================
// 🎓 TOUR GUIDE CONFIGURATION
// ==========================================
const tourData = {
    "Home": [
        { target: null, title: "Welcome Researcher!", text: "I am your interactive tour guide. Drug discovery is like finding a needle in a cosmic haystack. Let's learn how Quantum Computing speeds this up." },
        { target: "tour-home-grid", title: "The Sandbox", text: "We have built educational modules. 👉 <strong>Action:</strong> Click on the <strong>🏎️ Algorithm Race</strong> module in the grid." }
    ],
    "Algorithm Race": [
        { target: null, title: "Linear vs Exponential", text: "<strong>Why is discovery slow?</strong><br>Classical computers test one molecule at a time sequentially ($O(N)$). Quantum computers use superposition to test the entire database simultaneously ($O(\\sqrt{N})$)." },
        { target: "start-race-btn", title: "Run the Race", text: "👉 <strong>Action:</strong> Click <strong>Start Simulation Race</strong> to watch the visual bottleneck in real-time." }
    ],
    "Cost Simulator": [
        { target: null, title: "Economics of R&D", text: "<strong>Why is discovery expensive?</strong><br>Testing millions of molecules classically requires massive supercomputer time and funding." },
        { target: "cost-n", title: "Simulate Cost", text: "👉 <strong>Action:</strong> Change the <em>Number of compounds</em> box to <strong>10,000,000</strong> to see how standard pharmaceutical costs spiral out of control." }
    ],
    "Hilbert Space": [
        { target: "oracle-sidebar", title: "The Quantum Oracle", text: "This sidebar is our Quantum Oracle. We use it to map the physical properties a perfect drug should have." },
        { target: "tour-qed", title: "QED Score", text: "<strong>Quantitative Estimate of Druglikeness.</strong> A score closer to 1.0 means it looks like a real drug. 👉 <strong>Action:</strong> Try sliding it to 0.70." },
        { target: "tour-mw", title: "Molecular Weight", text: "<strong>Lipinski's Rule:</strong> Heavier molecules can't cross human cell walls. We cap it at 500 Daltons to ensure it's absorbable." },
        { target: "tour-logp", title: "LogP (Lipophilicity)", text: "This measures fat solubility. Too high, and it gets stuck in human fat tissue. Too low, and it cannot cross membranes." },
        { target: "tour-noise", title: "Decoherence Noise", text: "Real quantum hardware is noisy! Try adding 0.05 noise to see how error affects probability." },
        { target: "run-btn", title: "Execute Math", text: "👉 <strong>Action:</strong> Click <strong>Apply & Run</strong> to execute Grover's algorithm and isolate the targets!" }
    ],
    "Top Candidate": [
        { target: null, title: "Wave Collapse!", text: "The quantum probability collapsed into our best match! Look at the Physical Profile and 3D structure." },
        { target: "tc-lipinski-badges", title: "Lipinski Compliance", text: "Christopher Lipinski proved that oral drugs must obey strict rules. Our Oracle enforced this automatically!" }
    ],
    "Lipinski Filter": [
        { target: null, title: "The Rule of 5 Validation", text: "Here we batch-process our dataset against Lipinski's classical guardrails. Any molecule that fails more than 1 rule is flagged as Non-Compliant (toxic or unabsorbable)." }
    ],
    "Molecule Designer": [
        { target: "ketcher-frame", title: "Be the Chemist!", text: "Our Top Hit was injected into the canvas. Try drawing a new atom (like Oxygen 'O') onto the structure." },
        { target: "designer-smi", title: "Check Alignment", text: "👉 <strong>Action:</strong> Click <strong>Run Diagnostics</strong> and check if your new creation stays inside our active Oracle bounds on the left!" }
    ]
};

let currentTourStep = 0;
let currentTourSequence = [];

document.addEventListener("DOMContentLoaded", async () => {
    fetch3Dmol("CC1=C(C=C(C=C1)NC(=O)C2=CC=C(C=C2)CN3CCN(CC3)C)NC4=NC=CC(=N4)C5=CN=CC=C5", "home-3dmol");

    document.getElementById('theme-toggle').addEventListener('click', (e) => {
        const html = document.documentElement;
        if(html.getAttribute('data-theme') === 'light') { html.setAttribute('data-theme', 'dark'); e.target.innerText='☀️'; }
        else { html.setAttribute('data-theme', 'light'); e.target.innerText='🌙'; }
        if(globalSimData) renderSpecificView(currentActivePage);
    });

    await runQuantumSimulation();
});

function navigate(pageId) {
    currentActivePage = pageId;
    document.querySelectorAll('.nav-btn').forEach(btn => btn.classList.remove('active'));

    const caller = Array.from(document.querySelectorAll('.nav-btn')).find(b => b.innerText.includes(pageId) || b.getAttribute('onclick').includes(pageId));
    if(caller) caller.classList.add('active');

    document.querySelectorAll('.page-section').forEach(p => p.classList.remove('active'));
    document.querySelectorAll('.module-view').forEach(m => m.style.display = 'none');

    const sidebar = document.getElementById('oracle-sidebar');
    const appGrid = document.getElementById('main-app-grid');

    if (['Similarity Search', 'Lipinski Filter', 'Q-Visualizer', 'Algorithm Race', 'Cost Simulator'].includes(pageId)) {
        if(sidebar) sidebar.style.display = 'none';
        if(appGrid) appGrid.style.gridTemplateColumns = '1fr';
    } else {
        if(sidebar) sidebar.style.display = 'block';
        if(appGrid) appGrid.style.gridTemplateColumns = '1fr 3.5fr';
    }

    if(pageId === 'Home') {
        document.getElementById('Home').classList.add('active');
        const wrapper = document.getElementById('App-Wrapper');
        if (wrapper) wrapper.style.display = 'none';
    } else {
        const wrapper = document.getElementById('App-Wrapper');
        if (wrapper) { wrapper.style.display = 'block'; wrapper.classList.add('active'); }

        const titleEl = document.getElementById('current-module-title');
        if (titleEl) titleEl.innerText = pageId;

        const targetDiv = document.getElementById(pageId.replace(/\s+/g, '-'));
        if(targetDiv) targetDiv.style.display = 'block';
        setTimeout(() => { renderSpecificView(pageId); }, 150);
    }

    startTourSequence(pageId);
    window.scrollTo({ top: 0, behavior: 'smooth' });
}

async function runQuantumSimulation() {
    const btn = document.getElementById('run-btn');
    if(btn) btn.innerText = "⚛️ Calculating...";

    const statusHdr = document.getElementById('status-hdr');
    if(statusHdr) statusHdr.innerHTML = `<span style="color:var(--accent2-color);">Running Qiskit Math...</span>`;

    const config = {
        t_qed: parseFloat(document.getElementById('t_qed').value),
        min_mw: parseFloat(document.getElementById('min_mw').value), max_mw: parseFloat(document.getElementById('max_mw').value),
        min_logp: parseFloat(document.getElementById('min_logp').value), max_logp: parseFloat(document.getElementById('max_logp').value),
        min_tpsa: parseFloat(document.getElementById('min_tpsa').value), max_tpsa: parseFloat(document.getElementById('max_tpsa').value),
        max_hbd: parseInt(document.getElementById('max_hbd').value), max_hba: parseInt(document.getElementById('max_hba').value),
        max_rotbonds: parseInt(document.getElementById('max_rotbonds').value), max_rings: parseInt(document.getElementById('max_rings').value),
        grover_steps: parseInt(document.getElementById('grover_steps').value), noise_rate: parseFloat(document.getElementById('noise_rate').value)
    };

    try {
        const res = await fetch(`${API_URL}/simulate`, { method: 'POST', headers:{'Content-Type':'application/json'}, body: JSON.stringify(config) });
        if(!res.ok) throw new Error("Backend Crash");
        globalSimData = await res.json();

        lastLoadedSmiles = null;

        const hdrStr = `SEARCH SPACE: <span style="color:var(--text-color); font-weight:bold;">${globalSimData.N.toLocaleString()}</span> &nbsp;|&nbsp; ORACLE TARGETS: <span style="color:var(--accent-color); font-weight:bold;">${globalSimData.targets_count.toLocaleString()}</span> &nbsp;|&nbsp; ALLOCATED QUBITS: <span style="color:var(--accent-color); font-weight:bold;">${globalSimData.num_qubits}</span>`;
        if(statusHdr) statusHdr.innerHTML = hdrStr;

        const statN = document.getElementById('stat-N'); if(statN) statN.innerText = globalSimData.N.toLocaleString();
        const statDim = document.getElementById('stat-dim'); if(statDim) statDim.innerText = globalSimData.dim.toLocaleString();
        const statOpt = document.getElementById('stat-opt'); if(statOpt) statOpt.innerText = globalSimData.opt_steps;
        const hTarget = document.getElementById('home-target-count'); if(hTarget) hTarget.innerText = globalSimData.targets_count.toLocaleString();
        const aTarget = document.getElementById('app-target-count'); if(aTarget) aTarget.innerText = globalSimData.targets_count.toLocaleString();

        renderSpecificView(currentActivePage);
    } catch(err) {
        console.error(err);
        if(statusHdr) statusHdr.innerHTML = `<span style="color:#ef4444; font-weight:bold;">BACKEND ERROR</span>`;
    }

    if(btn) btn.innerText = "Apply & Run Simulation";
}

function getLayoutProps() {
    const isLight = document.documentElement.getAttribute('data-theme') === 'light';
    return { fontCol: isLight ? '#0f172a' : '#ffffff', bg: 'rgba(0,0,0,0)' };
}

function renderSpecificView(pageId) {
    if(!globalSimData || !globalSimData.display_df) return;
    if(pageId === 'Hilbert Space') renderHilbertSpace();
    if(pageId === 'Top Candidate') renderTopCandidate();
    if(pageId === 'Grover Analytics') renderAnalytics();
    if(pageId === 'Lipinski Filter') renderLipinskiFilter();
    if(pageId === 'Algorithm Race') initRaceView();
    if(pageId === 'Cost Simulator') renderCostSimulator(false);
    if(pageId === 'Molecule Designer') renderMoleculeDesigner();
}

function renderCostSimulator(animate = false) {
    const l = getLayoutProps();
    let N = parseInt(document.getElementById('cost-n').value) || 1000000;
    let costPer = parseFloat(document.getElementById('cost-per').value) || 10;
    let timePer = parseFloat(document.getElementById('cost-time').value) || 1;

    let cCost = N * costPer; let cTime = N * timePer;
    let qEst = Math.floor(Math.sqrt(N));
    let qCost = qEst * costPer; let qTime = qEst * timePer;
    let speedup = qTime > 0 ? cTime / qTime : 0;

    document.getElementById('metric-c-cost').innerText = `$${cCost.toLocaleString()}`;
    document.getElementById('metric-c-time').innerText = `${cTime.toLocaleString()} sec`;
    document.getElementById('metric-q-cost').innerText = `$${qCost.toLocaleString()}`;
    document.getElementById('metric-q-time').innerText = `${qTime.toLocaleString()} sec`;
    document.getElementById('metric-speedup').innerText = `${speedup.toLocaleString(undefined, {maximumFractionDigits: 0})}x`;
    document.getElementById('real-world-text').innerHTML = `Your classical simulation at N=${N.toLocaleString()} already costs <strong>$${cCost.toLocaleString()}</strong>, showing how early-stage screening contributes significantly to pharmaceutical R&D expenses. Quantum-inspired tools reduce this to $${qCost.toLocaleString()}.`;

    Plotly.newPlot('cost-bar-chart', [{ x: ["Classical", "Quantum-Inspired"], y: [cCost, qCost], type: 'bar', text: [`$${cCost.toLocaleString()}`, `$${qCost.toLocaleString()}`], textposition: 'auto', marker: { color: ['#ef4444', '#10b981'] } }], { title: {text: 'Total Cost Comparison ($)', font: {size: 14}}, paper_bgcolor: l.bg, plot_bgcolor: l.bg, font: { color: l.fontCol }, margin: {l: 50, r: 20, t: 40, b: 30} });

    let sizes = []; let maxN = Math.max(N, 1000);
    for(let i=1; i<=100; i++) { sizes.push(Math.floor((i/100) * maxN)); }
    let cCurve = sizes.map(s => s); let qCurve = sizes.map(s => Math.sqrt(s));

    Plotly.newPlot('cost-log-chart', [ { x: sizes, y: cCurve, mode: 'lines', name: "Classical O(N)", line: {color: '#ef4444', width: 3} }, { x: sizes, y: qCurve, mode: 'lines', name: "Quantum O(√N)", line: {color: '#10b981', width: 3} } ], { title: {text: 'Logarithmic View (Real Difference)', font: {size: 14}}, yaxis: {type: 'log'}, paper_bgcolor: l.bg, plot_bgcolor: l.bg, font: { color: l.fontCol }, margin: {l: 50, r: 20, t: 40, b: 30}, legend: {x: 0, y: 1} });

    if (!animate) { Plotly.newPlot('cost-scaling-chart', [ { x: sizes, y: cCurve, mode: 'lines', name: "Classical O(N)", line: {color: '#ef4444', width: 3} }, { x: sizes, y: qCurve, mode: 'lines', name: "Quantum O(√N)", line: {color: '#10b981', width: 3} } ], { paper_bgcolor: l.bg, plot_bgcolor: l.bg, font: { color: l.fontCol }, margin: {l: 50, r: 20, t: 20, b: 30}, legend: {x: 0, y: 1} }); }
}

function animateCostCurve() {
    if(costAnimationInterval) clearInterval(costAnimationInterval);
    let N = parseInt(document.getElementById('cost-n').value) || 1000000;
    let sizes = []; let maxN = Math.max(N, 1000);
    for(let i=1; i<=100; i++) { sizes.push(Math.floor((i/100) * maxN)); }
    let cCurve = sizes.map(s => s); let qCurve = sizes.map(s => Math.sqrt(s));

    let step = 5; const l = getLayoutProps();
    Plotly.newPlot('cost-scaling-chart', [ { x: sizes.slice(0, step), y: cCurve.slice(0, step), mode: 'lines', name: "Classical O(N)", line: {color: '#ef4444', width: 3} }, { x: sizes.slice(0, step), y: qCurve.slice(0, step), mode: 'lines', name: "Quantum O(√N)", line: {color: '#10b981', width: 3} } ], { xaxis: {range: [0, maxN]}, yaxis: {range: [0, maxN]}, paper_bgcolor: l.bg, plot_bgcolor: l.bg, font: { color: l.fontCol }, margin: {l: 50, r: 20, t: 20, b: 30}, legend: {x: 0, y: 1} });

    costAnimationInterval = setInterval(() => {
        step += 2;
        Plotly.update('cost-scaling-chart', { x: [sizes.slice(0, step), sizes.slice(0, step)], y: [cCurve.slice(0, step), qCurve.slice(0, step)] });
        if(step >= sizes.length) clearInterval(costAnimationInterval);
    }, 40);
}

function initRaceView() {
    document.getElementById('race-classic-text').innerText = `0 / ${globalSimData.N.toLocaleString()}`;
    document.getElementById('race-quantum-text').innerText = `0 / ${globalSimData.opt_steps}`;
    document.getElementById('race-classic-fill').style.width = '0%';
    document.getElementById('race-quantum-fill').style.width = '0%';
    document.getElementById('race-classic-cost').innerText = "$0.00";
    document.getElementById('race-quantum-cost').innerText = "$0.00";
    document.getElementById('race-winner-text').style.display = 'none';
}

function startRace() {
    if(!globalSimData) return;
    const btn = document.getElementById('start-race-btn');
    btn.disabled = true;
    initRaceView();
    const N = globalSimData.N; const steps = globalSimData.opt_steps;

    let qCount = 0;
    const qInterval = setInterval(() => {
        qCount++; const pct = (qCount / steps) * 100;
        document.getElementById('race-quantum-fill').style.width = `${pct}%`;
        document.getElementById('race-quantum-text').innerText = `${qCount} / ${steps}`;
        document.getElementById('race-quantum-cost').innerText = "$" + (qCount * 0.05).toFixed(2);
        if(qCount >= steps) { clearInterval(qInterval); document.getElementById('race-winner-text').style.display = 'block'; }
    }, 1000 / steps);

    let cCount = 0;
    const cInterval = setInterval(() => {
        cCount += Math.floor(N / 80); if(cCount > N) cCount = N;
        const pct = (cCount / N) * 100;
        document.getElementById('race-classic-fill').style.width = `${pct}%`;
        document.getElementById('race-classic-text').innerText = `${cCount.toLocaleString()} / ${N.toLocaleString()}`;
        document.getElementById('race-classic-cost').innerText = "$" + (cCount * 0.01).toFixed(2);
        if(cCount >= N) { clearInterval(cInterval); btn.disabled = false; }
    }, 100);
}

function renderHilbertSpace() {
    if(document.getElementById('Hilbert-Space').style.display === 'none') return;
    const data = globalSimData.display_df;
    const l = getLayoutProps();
    const probValues = data.map(d => d.Probability);
    let minP = Math.min(...probValues); let maxP = Math.max(...probValues);
    let cminVal = minP; let cmaxVal = maxP;
    if ((maxP - minP) < 1e-5) { cminVal = 0; cmaxVal = maxP > 0 ? maxP * 3 : 1; }

    const sizes = data.map(d => d.Classification === '★ MEASURED WAVE COLLAPSE' ? 14 : 5);
    const symbols = data.map(d => d.Classification === '★ MEASURED WAVE COLLAPSE' ? 'diamond' : 'circle');
    const lineColors = data.map(d => d.Classification === '★ MEASURED WAVE COLLAPSE' ? '#ffffff' : 'transparent');
    const lineWidths = data.map(d => d.Classification === '★ MEASURED WAVE COLLAPSE' ? 2 : 0);

    Plotly.newPlot('hilbert-plot', [{
        x: data.map(d=>d.MW), y: data.map(d=>d.LogP), z: probValues, text: data.map(d=>d.Name),
        mode: 'markers', type: 'scatter3d',
        marker: { size: sizes, color: probValues, cmin: cminVal, cmax: cmaxVal, colorscale: [[0, '#1e293b'], [0.1, '#3b82f6'], [0.5, '#a855f7'], [1.0, '#f59e0b']], showscale: true, colorbar: { title: '|Ψ|² Prob', thickness: 15, tickfont: {color: l.fontCol}, titlefont: {color: l.fontCol, size: 12} }, opacity: 0.9, symbol: symbols, line: { color: lineColors, width: lineWidths } }
    }], { paper_bgcolor: l.bg, plot_bgcolor: l.bg, font: { color: l.fontCol }, scene: { xaxis: {title: 'Mass (Da)'}, yaxis: {title: 'LogP'}, zaxis: { title: '|Ψ|² Probability', range: [0, maxP > 0 ? maxP * 1.5 : 1] } }, margin: {l:0,r:0,t:0,b:0} });

    let tableHTML = `<table class="custom-table"><tr><th>Name</th><th>Classification</th><th>MW</th><th>LogP</th><th>QED</th><th>TPSA</th><th>HBD</th><th>HBA</th><th>Probability</th></tr>`;
    data.slice(0, 20).forEach(d => { tableHTML += `<tr><td>${d.Name}</td><td style="color:${d.Classification==='★ MEASURED WAVE COLLAPSE'?'#f59e0b':'#3b82f6'}">${d.Classification}</td><td>${d.MW.toFixed(2)}</td><td>${d.LogP.toFixed(2)}</td><td>${d.QED.toFixed(3)}</td><td>${d.TPSA.toFixed(1)}</td><td>${d.HBD}</td><td>${d.HBA}</td><td>${(d.Probability*100).toFixed(4)}%</td></tr>`; });
    tableHTML += `</table>`;
    const tbContainer = document.getElementById('hilbert-table-container');
    if (tbContainer) tbContainer.innerHTML = tableHTML;
}

function renderTopCandidate() {
    if(document.getElementById('Top-Candidate').style.display === 'none') return;
    const hit = globalSimData.top_hit;
    document.getElementById('tc-name').innerText = hit.Name;
    document.getElementById('tc-prob').innerText = (hit.Probability * 100).toFixed(4) + "%";
    document.getElementById('tc-mw').innerText = hit.MW.toFixed(2) + " Da";
    document.getElementById('tc-logp').innerText = hit.LogP.toFixed(2);
    document.getElementById('tc-qed').innerText = hit.QED.toFixed(3);
    document.getElementById('tc-tpsa').innerText = hit.TPSA.toFixed(1) + " Å²";

    let v = 0; let badgeHtml = '';
    const checks = [{name: "MW ≤ 500 Da", pass: hit.MW <= 500}, {name: "LogP ≤ 5", pass: hit.LogP <= 5}, {name: "HBD ≤ 5", pass: hit.HBD <= 5}, {name: "HBA ≤ 10", pass: hit.HBA <= 10}];
    checks.forEach(c => { if(!c.pass) v++; badgeHtml += `<div class="${c.pass ? 'badge-pass' : 'badge-fail'}" style="margin-bottom: 8px;">${c.pass ? 'PASS' : 'FAIL'} &mdash; ${c.name}</div>`; });
    document.getElementById('tc-lipinski-badges').innerHTML = badgeHtml;
    document.getElementById('tc-verdict').innerHTML = `<div class="card-header" style="font-size: 0.75rem; border:none; margin-bottom:5px;">SYSTEM VERDICT:</div><div class="${v<=1?'badge-pass':'badge-fail'}">${v<=1?'COMPLIANT':`NON-COMPLIANT (${v})`}</div>`;

    const l = getLayoutProps();
    Plotly.newPlot('radar-plot', [{ type: 'scatterpolar', r: [hit.QED, Math.min(hit.MW/500,1), Math.min((hit.LogP+5)/15,1), Math.min(hit.TPSA/140,1), Math.min(hit.HBD/5,1), Math.min(hit.HBA/10,1)], theta: ['QED','MW','LogP','TPSA','HBD','HBA'], fill: 'toself', line: {color: '#3b82f6'} }], { polar: {radialaxis:{visible:false}, bgcolor: l.bg}, paper_bgcolor: l.bg, font:{color:l.fontCol}, margin:{l:30,r:30,t:10,b:10} });
    if(hit.SMILES && hit.SMILES !== "C") fetch3Dmol(hit.SMILES, 'tc-3dmol');
}

function renderLipinskiFilter() {
    if(document.getElementById('Lipinski-Filter').style.display === 'none') return;
    const l = getLayoutProps();
    let passCount = 0; let failCount = 0;
    const realDrugs = globalSimData.display_df.filter(d => d.Is_Real).slice(0, 100);

    let tableHTML = `<table class="custom-table"><tr><th>Name</th><th>MW</th><th>LogP</th><th>HBD</th><th>HBA</th><th>Violations</th><th>Status</th></tr>`;
    realDrugs.forEach(d => {
        let v = 0; if(d.MW > 500) v++; if(d.LogP > 5) v++; if(d.HBD > 5) v++; if(d.HBA > 10) v++;
        const status = v <= 1 ? "PASS" : "FAIL";
        if(status === "PASS") passCount++; else failCount++;
        tableHTML += `<tr><td>${d.Name}</td><td>${d.MW.toFixed(2)}</td><td>${d.LogP.toFixed(2)}</td><td>${d.HBD}</td><td>${d.HBA}</td><td>${v}</td><td><span class="${status==='PASS'?'badge-pass':'badge-fail'}">${status}</span></td></tr>`;
    });
    tableHTML += `</table>`;
    document.getElementById('lipinski-table-container').innerHTML = tableHTML;
    Plotly.newPlot('lipinski-pie', [{ values: [passCount, failCount], labels: ['Compliant', 'Non-Compliant'], type: 'pie', hole: 0.6, marker: { colors: ['#3b82f6', '#f59e0b'] } }], { title: 'Dataset Compliance', paper_bgcolor: l.bg, font: { color: l.fontCol }, margin: {l:20, r:20, t:40, b:20} });
}

function renderAnalytics() {
    if(document.getElementById('Grover-Analytics').style.display === 'none') return;
    const l = getLayoutProps();
    const trace = globalSimData.psi_history.map(stepProbs => Math.max(...stepProbs));
    const xSteps = Array.from({length: trace.length}, (_, i) => i);

    Plotly.newPlot('grover-line', [{ x: xSteps, y: trace, type: 'scatter', mode: 'lines+markers', line:{color:'#3b82f6', width: 3}, marker:{size:8} }], { title: {text: 'Amplitude Amplification Trace (Qiskit)', font: {size: 14}}, paper_bgcolor: l.bg, plot_bgcolor: l.bg, font:{color:l.fontCol}, xaxis:{title:'Grover Iteration', tickmode:'linear'}, yaxis:{title:'|Ψ|² Probability'}, margin:{l:50,r:20,t:40,b:40} });
    const top50 = globalSimData.display_df.slice(0,50);
    Plotly.newPlot('grover-bar', [{ x: top50.map((_,i)=>i), y: top50.map(d=>d.Probability), type: 'bar', marker:{color: top50.map(d=>d.Classification==='Suppressed Noise'? '#475569':'#3b82f6')} }], { title: {text:'Top 50 Qubit States', font: {size: 14}}, xaxis:{title:'State Rank'}, yaxis:{title:'|Ψ|² Probability'}, paper_bgcolor: l.bg, plot_bgcolor: l.bg, font:{color:l.fontCol}, margin:{l:50,r:20,t:40,b:40} });

    document.getElementById('log-init').innerHTML = `Allocated ${globalSimData.num_qubits} Qubits dynamically.<br>Hilbert space dimension 2^${globalSimData.num_qubits} = ${globalSimData.dim.toLocaleString()}.`;
    document.getElementById('log-oracle').innerHTML = `Phase-flip applied. Target count: ${globalSimData.targets_count.toLocaleString()}`;
    const nr = document.getElementById('noise_rate').value;
    document.getElementById('log-noise').innerHTML = `Noise depolarizing rate: ${nr}.<br>Coherence preserved: ${((1-nr)*100).toFixed(0)}%`;
    document.getElementById('log-speedup').innerHTML = `O(N) = ${globalSimData.N.toLocaleString()} iters vs O(√N) = ${globalSimData.opt_steps} iters.`;
    document.getElementById('log-collapse').innerHTML = `Collapsed state:<br><b>${globalSimData.top_hit.Name}</b> at ${(globalSimData.top_hit.Probability*100).toFixed(4)}% amp.`;
}

async function runSimilarity() {
    const smi = document.getElementById('sim-input').value;
    if(!smi) return;
    const res = await fetch(`${API_URL}/similarity`, { method: 'POST', headers:{'Content-Type':'application/json'}, body: JSON.stringify({smiles: smi}) });
    const data = await res.json();
    const l = getLayoutProps();
    Plotly.newPlot('sim-plot', [{ x: data.matches.map(d=>d.Name), y: data.matches.map(d=>d.Tanimoto), type: 'bar', marker:{color:'#3b82f6'} }], { paper_bgcolor: l.bg, plot_bgcolor: l.bg, font:{color:l.fontCol}, margin:{l:40,r:20,t:20,b:40} });
    document.getElementById('sim-top-match').innerText = `Top Match: ${data.matches[0].Name} (${data.matches[0].Tanimoto.toFixed(3)})`;
    fetch3Dmol(data.matches[0].SMILES, 'sim-3dmol');
}

function renderMoleculeDesigner() {
    if(document.getElementById('Molecule-Designer').style.display === 'none') return;
    if(!globalSimData || !globalSimData.top_hit || !globalSimData.top_hit.SMILES) return;
    const targetSmiles = globalSimData.top_hit.SMILES;

    if (lastLoadedSmiles !== targetSmiles && targetSmiles !== "C") {
        const iframe = document.getElementById('ketcher-frame');
        if (iframe && iframe.contentWindow && iframe.contentWindow.ketcher) {
            iframe.contentWindow.ketcher.setMolecule(targetSmiles).then(() => {
                lastLoadedSmiles = targetSmiles;
                document.getElementById('designer-smi').value = targetSmiles;
                runDesignerDiagnostics();
            }).catch(e => console.log("Ketcher initialization waiting..."));
        } else {
            setTimeout(renderMoleculeDesigner, 500);
        }
    }
}

async function runDesignerDiagnostics() {
    let smi = document.getElementById('designer-smi').value;
    const iframe = document.getElementById('ketcher-frame');
    if (iframe && iframe.contentWindow && iframe.contentWindow.ketcher) {
        try {
            const kSmiles = await iframe.contentWindow.ketcher.getSmiles();
            if (kSmiles) { smi = kSmiles; document.getElementById('designer-smi').value = smi; }
        } catch(e) { console.log("Fetching Ketcher SMILES..."); }
    }
    if(!smi) return;
    try {
        const res = await fetch(`${API_URL}/molecule-props`, { method: 'POST', headers:{'Content-Type':'application/json'}, body: JSON.stringify({smiles: smi}) });
        const data = await res.json();
        const b = { qed: parseFloat(document.getElementById('t_qed').value), min_mw: parseFloat(document.getElementById('min_mw').value), max_mw: parseFloat(document.getElementById('max_mw').value), min_log: parseFloat(document.getElementById('min_logp').value), max_log: parseFloat(document.getElementById('max_logp').value), min_tpsa: parseFloat(document.getElementById('min_tpsa').value), max_tpsa: parseFloat(document.getElementById('max_tpsa').value), hbd: parseInt(document.getElementById('max_hbd').value), hba: parseInt(document.getElementById('max_hba').value), rot: parseInt(document.getElementById('max_rotbonds').value), ring: parseInt(document.getElementById('max_rings').value) };
        const inOracle = (data.props.QED >= b.qed && data.props.MW >= b.min_mw && data.props.MW <= b.max_mw && data.props.LogP >= b.min_log && data.props.LogP <= b.max_log && data.props.TPSA >= b.min_tpsa && data.props.TPSA <= b.max_tpsa && data.props.HBD <= b.hbd && data.props.HBA <= b.hba && data.props.RotBonds <= b.rot && data.props.Rings <= b.ring);

        let html = `
            <div style='display: flex; justify-content: space-between; margin-bottom: 0.5rem;'><div class='data-lbl'>Mol Weight</div><div style='font-family: JetBrains Mono; font-weight: 700; font-size: 1.1rem;'>${data.props.MW.toFixed(2)} Da</div></div>
            <div style='display: flex; justify-content: space-between; margin-bottom: 0.5rem;'><div class='data-lbl'>LogP</div><div style='font-family: JetBrains Mono; font-weight: 700; font-size: 1.1rem;'>${data.props.LogP.toFixed(2)}</div></div>
            <div style='display: flex; justify-content: space-between; margin-bottom: 0.5rem;'><div class='data-lbl'>TPSA</div><div style='font-family: JetBrains Mono; font-weight: 700; font-size: 1.1rem;'>${data.props.TPSA.toFixed(1)} Å²</div></div>
            <div style='display: flex; justify-content: space-between; margin-bottom: 0.5rem;'><div class='data-lbl'>QED Score</div><div style='font-family: JetBrains Mono; font-weight: 700; font-size: 1.1rem;'>${data.props.QED.toFixed(3)}</div></div>
            <div style='display: flex; justify-content: space-between; margin-bottom: 0.5rem;'><div class='data-lbl'>HB Donors</div><div style='font-family: JetBrains Mono; font-weight: 700; font-size: 1.1rem;'>${data.props.HBD}</div></div>
            <div style='display: flex; justify-content: space-between; margin-bottom: 1.5rem;'><div class='data-lbl'>HB Acceptors</div><div style='font-family: JetBrains Mono; font-weight: 700; font-size: 1.1rem;'>${data.props.HBA}</div></div>
            <div class='card-header' style='font-size: 0.75rem; border:none; margin-bottom: 5px;'>LIPINSKI RULE OF 5</div>
        `;
        for (const [rule, pass] of Object.entries(data.checks)) { html += `<div class="${pass ? 'badge-pass' : 'badge-fail'}" style="margin-bottom: 5px;">${pass?'PASS':'FAIL'} - ${rule}</div>`; }
        html += `<div class='card-header' style='font-size: 0.75rem; border:none; margin-top: 1.5rem; margin-bottom: 5px;'>ORACLE ALIGNMENT</div>`;
        html += `<div class="${inOracle ? 'badge-pass' : 'badge-fail'}">${inOracle?'✓ IN ACTIVE BOUNDS':'✗ OUTSIDE ACTIVE BOUNDS'}</div>`;
        document.getElementById('designer-results').innerHTML = html;
        fetch3Dmol(smi, 'designer-3dmol');
    } catch { document.getElementById('designer-results').innerHTML = `<div class="badge-fail">Invalid SMILES Structure</div>`; }
}

async function runQVisualizer() {
    const statBox = document.getElementById('qv-status-box');
    const stat = document.getElementById('qv-status');
    statBox.style.display = 'block';
    const N = 64; let amps = Array(N).fill(1/Math.sqrt(N));
    const targetIdx = 32; const l = getLayoutProps();
    const xLabels = Array.from({length: N}, (_, i) => i === targetIdx ? "ASPIRIN" : `SYN-${i.toString().padStart(3, '0')}`);

    function plotIt(ydata, title, yTitle="Probability Amplitude", yRange=[-0.4, 0.8]) {
        Plotly.newPlot('qv-plot', [{ x: xLabels, y: ydata, type: 'bar', marker:{color: ydata.map((_,i)=> i===targetIdx ? '#f59e0b' : '#3b82f6')} }], { title: title, paper_bgcolor: l.bg, plot_bgcolor: l.bg, font:{color:l.fontCol}, yaxis:{title: yTitle, range:yRange}, xaxis: {showticklabels: false}, margin:{l:50,r:20,t:40,b:40} });
    }

    stat.innerHTML = "<b>PHASE 1: UNIFORM SUPERPOSITION</b><br>The Statevector is initialized. Hadamard gates apply a uniform transformation across the basis states. Every molecule exists simultaneously with an identical initial amplitude.";
    plotIt(amps, "Step 1: Uniform Superposition", "Probability Amplitude", [-0.4, 0.4]);
    await new Promise(r => setTimeout(r, 2000));

    const steps = Math.floor((Math.PI / 4) * Math.sqrt(N));
    for(let step=0; step<steps; step++) {
        stat.innerHTML = `<b>ITERATION ${step+1}: ORACLE PHASE FLIP</b><br>The unitary Diagonal Oracle matrix maps the constraints and applies a negative phase shift exclusively to the target index (Aspirin).`;
        amps[targetIdx] *= -1;
        plotIt(amps, `Iteration ${step+1}: Oracle Phase Flip`, "Probability Amplitude", [-0.5, 1.0]);
        await new Promise(r => setTimeout(r, 1500));

        stat.innerHTML = `<b>ITERATION ${step+1}: GROVER DIFFUSION</b><br>The Diffuser applies inversion-about-the-mean. Notice how the target amplitude surges while the noise amplitudes are suppressed via destructive interference.`;
        const mean = amps.reduce((a,b)=>a+b,0)/N;
        amps = amps.map(a => 2*mean - a);
        plotIt(amps, `Iteration ${step+1}: Amplitude Amplification`, "Probability Amplitude", [-0.5, 1.0]);
        await new Promise(r => setTimeout(r, 1500));
    }

    stat.innerHTML = "<b>PHASE 3: STATE MEASUREMENT (WAVE COLLAPSE)</b><br>Probability extraction completes. The vector collapses into a deterministic classical output. Aspirin is isolated perfectly.";
    plotIt(amps.map(a=>a*a), "Final State: Wave Collapse Measurement", "|Ψ|² Probability", [0, 1.05]);
}

async function fetch3Dmol(smiles, elementId) {
    try {
        const res = await fetch(`${API_URL}/3dmol?smiles=${encodeURIComponent(smiles)}`);
        const data = await res.json();
        const container = document.getElementById(elementId);
        container.innerHTML = "";
        let viewer = $3Dmol.createViewer(elementId, {backgroundColor: 'rgba(0,0,0,0)'});
        viewer.addModel(data.mol_block, "sdf");
        viewer.setStyle({}, {stick:{colorscheme:'blueWhiteCarbon', radius:0.2}});
        viewer.zoomTo(); viewer.zoom(1.8); viewer.spin(true); viewer.render();
    } catch(err) { console.error("3Dmol Error"); }
}

// ==========================================
// 🎓 STUDENT / EARLY RESEARCHER TOUR GUIDE
// ==========================================
function toggleStudentMode() {
    studentModeActive = !studentModeActive;
    const btn = document.getElementById('student-mode-toggle');
    const assistant = document.getElementById('student-assistant');

    if(studentModeActive) {
        btn.innerText = '🎓 Mode: Student';
        btn.style.background = 'var(--accent2-color)';
        btn.style.color = '#fff';
        assistant.style.display = 'block';
        assistant.animate([{opacity: 0, transform: 'translateY(20px)'}, {opacity: 1, transform: 'translateY(0)'}], {duration: 300, fill: 'forwards'});
        startTourSequence(currentActivePage);
    } else {
        btn.innerText = '🎓 Mode: Experienced';
        btn.style.background = 'transparent';
        btn.style.color = 'var(--accent2-color)';
        assistant.style.display = 'none';
        clearAllHighlights();
    }
}

function startTourSequence(pageId) {
    if(!studentModeActive) return;

    // Default fallback if page isn't in tourData
    currentTourSequence = tourData[pageId] || [{target: null, title: "Explore freely", text: "Use the navigation above to select a module."}];
    currentTourStep = 0;
    renderTourStep();
}

function renderTourStep() {
    if(!studentModeActive) return;
    clearAllHighlights();

    const stepInfo = currentTourSequence[currentTourStep];
    const textEl = document.getElementById('assistant-text');
    const titleEl = document.getElementById('tour-title');
    const progressEl = document.getElementById('tour-progress');

    // Add glowing highlight if a target is defined
    if(stepInfo.target) {
        const targetEl = document.getElementById(stepInfo.target);
        if(targetEl) targetEl.classList.add('student-focus');
    }

    textEl.style.opacity = 0;
    setTimeout(() => {
        titleEl.innerText = stepInfo.title;
        textEl.innerHTML = stepInfo.text;
        progressEl.innerText = `${currentTourStep + 1}/${currentTourSequence.length}`;
        textEl.style.opacity = 1;
        textEl.style.transition = "opacity 0.3s";
    }, 150);
}

function tourNext() {
    if (currentTourStep < currentTourSequence.length - 1) {
        currentTourStep++;
        renderTourStep();
    }
}

function tourPrev() {
    if (currentTourStep > 0) {
        currentTourStep--;
        renderTourStep();
    }
}

function clearAllHighlights() {
    // Remove glow from all possible targets
    document.querySelectorAll('.student-focus').forEach(el => el.classList.remove('student-focus'));
}