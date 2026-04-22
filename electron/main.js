const { app, BrowserWindow, ipcMain, dialog } = require("electron");
const path = require("path");
const fs = require("fs");
const os = require("os");
const { spawn, spawnSync } = require("child_process");

// -----------------------------------------------------------------------
// Binary resolution: WSL on Windows, system PATH otherwise
// -----------------------------------------------------------------------

const IS_WIN = process.platform === "win32";

// Resolve a bundled binary or fall back to system PATH
function resolveBin(name) {
    const platform = process.platform === "win32" ? "win"
        : process.platform === "darwin" ? "mac" : "linux";
    const local = path.join(__dirname, "bin", platform, name);
    if (fs.existsSync(local)) {
        try { fs.chmodSync(local, 0o755); } catch { }
        return local;
    }
    return name;
}

// -----------------------------------------------------------------------
// App window
// -----------------------------------------------------------------------
function createWindow() {
    const win = new BrowserWindow({
        width: 1280,
        height: 900,
        minWidth: 800,
        minHeight: 600,
        title: "SPLACE",
        webPreferences: {
            preload: path.join(__dirname, "preload.js"),
            contextIsolation: true,
            nodeIntegration: false,
        },
    });
    win.loadFile(path.join(__dirname, "..", "docs", "index.html"));
}

app.whenReady().then(() => {
    createWindow();
    app.on("activate", () => {
        if (BrowserWindow.getAllWindows().length === 0) createWindow();
    });
});

app.on("window-all-closed", () => {
    if (process.platform !== "darwin") app.quit();
});

// -----------------------------------------------------------------------
// IPC: CPU count
// -----------------------------------------------------------------------
ipcMain.handle("get-cpu-count", () => os.cpus().length);

// -----------------------------------------------------------------------
// IPC: run MAFFT alignment
// { files: {geneName: fastaContent}, params: string[], threads: number }
// -----------------------------------------------------------------------
ipcMain.on("run-analysis", async (event, { files, params, threads }) => {
    const sender = event.sender;
    const send = (channel, data) => { if (!sender.isDestroyed()) sender.send(channel, data); };

    const outputDir = path.join(os.homedir(), "Documents", "SPLACE", "sequences");
    try { fs.rmSync(outputDir, { recursive: true, force: true }); } catch {}
    fs.mkdirSync(outputDir, { recursive: true });

    const tmpDir = fs.mkdtempSync(path.join(os.tmpdir(), "splace-"));
    const markers = Object.keys(files);
    let done = 0, aligned = 0;
    const markerResults = new Array(markers.length).fill(null);

    // Parallel pool — at most `threads` MAFFT jobs running simultaneously
    const concurrency = Math.max(1, threads);
    let nextIdx = 0;

    async function worker() {
        while (true) {
            const idx = nextIdx++;
            if (idx >= markers.length) break;
            const name = markers[idx];
            const inputFile = path.join(tmpDir, `${name}.fasta`);
            const outputFile = path.join(outputDir, `${name}_aligned.fasta`);
            const rawContent = files[name];

            fs.writeFileSync(inputFile, rawContent, "utf8");

            send("analysis-progress", {
                marker: name, status: "running",
                message: `[${name}] Running MAFFT…`,
                done, total: markers.length,
            });

            const rawSeqs = parseFastaStats(rawContent);
            const { success, alignedContent, error } = await runMafft(
                inputFile, outputFile, params, threads,
                (msg) => send("analysis-progress", { message: `[${name}] ${msg}`, done, total: markers.length })
            );

            done++;
            if (success) aligned++;
            const alignedStats = success ? parseFastaStats(alignedContent) : null;

            send("analysis-progress", {
                marker: name,
                status: success ? "done" : "error",
                info: success ? `${alignedStats.numSeqs} seq · ${alignedStats.length} bp` : error,
                message: success ? `[${name}] ✓ Saved → ${outputFile}` : `[${name}] ✗ ${error}`,
                done, total: markers.length,
            });

            markerResults[idx] = {
                marker: name, rawContent,
                alignedContent: success ? alignedContent : null,
                rawStats: rawSeqs, alignedStats,
                trimmedContent: null, trimmedStats: null,
                outputFile: success ? outputFile : null,
            };
        }
    }

    await Promise.all(Array.from({ length: concurrency }, () => worker()));
    try { fs.rmSync(tmpDir, { recursive: true, force: true }); } catch {}

    const finalResults = markerResults.filter(Boolean);
    send("analysis-done", {
        phase: "mafft",
        success: aligned > 0,
        aligned, total: markers.length,
        outputDir,
        markerResults: finalResults,
    });
    _lastMarkerResults = finalResults;
});

// -----------------------------------------------------------------------
// IPC: run trimAl
// { markers: string[], params: string[] }
// (uses already-aligned files from previous run stored in markerResults)
// -----------------------------------------------------------------------
let _lastMarkerResults = [];

ipcMain.on("run-trimal", async (event, { markers, params }) => {
    const sender = event.sender;
    const send = (channel, data) => { if (!sender.isDestroyed()) sender.send(channel, data); };

    let done = 0, trimmed = 0;
    const updatedResults = _lastMarkerResults.map(r => ({ ...r }));
    const concurrency = Math.max(1, os.cpus().length - 2);
    let nextIdx = 0;

    async function worker() {
        while (true) {
            const idx = nextIdx++;
            if (idx >= markers.length) break;
            const name = markers[idx];
            const mr = updatedResults.find(r => r.marker === name);
            if (!mr || !mr.alignedContent || !mr.outputFile) {
                done++;
                send("analysis-progress", { marker: name, status: "error", message: `[${name}] No aligned file found`, done, total: markers.length });
                continue;
            }

            const inputFile  = mr.outputFile;
            const trimmedFile = inputFile.replace(/_aligned\.fasta$/, "_trimmed.fasta");

            send("analysis-progress", {
                marker: name, status: "running",
                message: `[${name}] Running trimAl…`,
                done, total: markers.length,
            });

            const { success, trimmedContent, error } = await runTrimal(
                inputFile, trimmedFile, params,
                (msg) => send("analysis-progress", { message: `[${name}] ${msg}`, done, total: markers.length })
            );

            done++;
            if (success) {
                trimmed++;
                mr.trimmedContent = trimmedContent;
                mr.trimmedStats = parseFastaStats(trimmedContent);
            }

            send("analysis-progress", {
                marker: name,
                status: success ? "done" : "error",
                info: success ? `${mr.trimmedStats.numSeqs} seq · ${mr.trimmedStats.length} bp` : error,
                message: success ? `[${name}] ✓ Trimmed → ${trimmedFile}` : `[${name}] ✗ ${error}`,
                done, total: markers.length,
            });
        }
    }

    await Promise.all(Array.from({ length: concurrency }, () => worker()));

    send("analysis-done", {
        phase: "trimal",
        success: trimmed > 0,
        aligned: trimmed, total: markers.length,
        markerResults: updatedResults,
    });
});

// -----------------------------------------------------------------------
// Helper: run a single MAFFT alignment
// On Windows: pipe FASTA via stdin to `wsl bash -c "mafft ... - 2>/dev/null"`
//   • avoids Windows paths entirely inside WSL
//   • avoids /dev/stderr WSL interop crash (exec 3>/dev/stderr in mafft script)
// On Linux/Mac: spawn mafft directly with a file path
// -----------------------------------------------------------------------
function runMafft(inputFile, outputFile, params, threads, onLog) {
    return new Promise((resolve) => {
        const args = [...params, "--thread", String(threads), "--quiet"];
        let stdout = "", proc;

        if (IS_WIN) {
            const rawContent = fs.readFileSync(inputFile, "utf8");
            const tmpId = Date.now();
            const wslIn = `/tmp/splace_in_${tmpId}.fasta`;
            // Write FASTA to WSL /tmp via stdin, then run MAFFT on the native WSL path.
            // stdio 'ignore' for stderr → fd 2 becomes /dev/null inside WSL,
            // so MAFFT's  exec 3>/dev/stderr  succeeds and the script runs normally.
            const shellCmd =
                `cat > '${wslIn}' && ` +
                `mafft ${args.map(a => `'${a}'`).join(" ")} '${wslIn}'; ` +
                `rm -f '${wslIn}'`;
            proc = spawn("wsl", ["bash", "-c", shellCmd], {
                shell: false,
                stdio: ["pipe", "pipe", "ignore"],
            });
            proc.stdin.write(rawContent, "utf8");
            proc.stdin.end();
        } else {
            proc = spawn(resolveBin("mafft"), [...args, inputFile], { shell: true });
            proc.stderr.on("data", (d) => {
                const msg = d.toString().trim();
                if (msg) onLog(msg);
            });
        }

        proc.stdout.on("data", (d) => { stdout += d.toString(); });

        proc.on("close", (code) => {
            if (code === 0 && stdout.trim()) {
                try {
                    fs.writeFileSync(outputFile, stdout, "utf8");
                    resolve({ success: true, alignedContent: stdout });
                } catch (e) {
                    resolve({ success: false, error: `Write failed: ${e.message}` });
                }
            } else {
                resolve({ success: false, error: `Exit code ${code}` });
            }
        });

        proc.on("error", (err) => {
            resolve({ success: false, error: err.message });
        });
    });
}

// -----------------------------------------------------------------------
// Helper: run a single trimAl trimming
// On Windows: write aligned FASTA via stdin to WSL /tmp, run trimal there,
//   read output back via cat — avoids /dev/stderr interop crash
// -----------------------------------------------------------------------
function runTrimal(inputFile, outputFile, params, onLog) {
    return new Promise((resolve) => {
        let stdout = "", proc;

        if (IS_WIN) {
            const rawContent = fs.readFileSync(inputFile, "utf8");
            const tmpId = Date.now();
            const wslIn = `/tmp/splace_in_${tmpId}.fasta`;
            const wslOut = `/tmp/splace_out_${tmpId}.fasta`;
            const paramsStr = params.map(a => `'${a}'`).join(" ");
            const shellCmd =
                `cat > '${wslIn}' && ` +
                `trimal -in '${wslIn}' -out '${wslOut}' ${paramsStr} && ` +
                `cat '${wslOut}'; rm -f '${wslIn}' '${wslOut}'`;
            proc = spawn("wsl", ["bash", "-c", shellCmd], {
                shell: false,
                stdio: ["pipe", "pipe", "ignore"],
            });
            proc.stdin.write(rawContent, "utf8");
            proc.stdin.end();
            proc.stdout.on("data", (d) => { stdout += d.toString(); });
            proc.on("close", (code) => {
                if (code === 0 && stdout.trim()) {
                    fs.writeFileSync(outputFile, stdout, "utf8");
                    resolve({ success: true, trimmedContent: stdout });
                } else {
                    resolve({ success: false, error: `Exit code ${code}` });
                }
            });
        } else {
            const args = ["-in", inputFile, "-out", outputFile, ...params];
            proc = spawn(resolveBin("trimal"), args, { shell: false });
            proc.stderr.on("data", (d) => {
                const msg = d.toString().trim();
                if (msg) onLog(msg);
            });
            proc.on("close", (code) => {
                if (code === 0 && fs.existsSync(outputFile)) {
                    const trimmedContent = fs.readFileSync(outputFile, "utf8");
                    resolve({ success: true, trimmedContent });
                } else {
                    resolve({ success: false, error: `Exit code ${code}` });
                }
            });
        }

        proc.on("error", (err) => {
            resolve({ success: false, error: err.message });
        });
    });
}


// -----------------------------------------------------------------------
// Helper: parse basic FASTA stats
// -----------------------------------------------------------------------
function parseFastaStats(text) {
    const seqs = [];
    let cur = null;
    for (const line of (text || "").split(/\r?\n/)) {
        if (line.startsWith(">")) {
            if (cur !== null) seqs.push(cur);
            cur = 0;
        } else if (cur !== null) {
            cur += line.trim().length;
        }
    }
    if (cur !== null) seqs.push(cur);

    const numSeqs = seqs.length;
    const length = seqs[0] || 0;
    const avgLen = numSeqs > 0 ? Math.round(seqs.reduce((a, b) => a + b, 0) / numSeqs) : 0;

    // Gap percentage of first (aligned) sequence set
    const gapCount = (text || "").replace(/^>.*$/gm, "").replace(/\n/g, "").split("").filter(c => c === "-").length;
    const totalChars = (text || "").replace(/^>.*$/gm, "").replace(/\n/g, "").length;
    const gapPct = totalChars > 0 ? ((gapCount / totalChars) * 100).toFixed(1) : "0.0";

    return { numSeqs, length, avgLen, gapPct };
}

// -----------------------------------------------------------------------
// IPC: run IQ-TREE3
// { nexus, partition, params, threads, outgroup }
// -----------------------------------------------------------------------
ipcMain.on("run-iqtree", async (event, { nexus, partition, params, threads, outgroup, perGene, geneFiles }) => {
    const sender = event.sender;
    const send = (channel, data) => { if (!sender.isDestroyed()) sender.send(channel, data); };

    const outputDir = path.join(os.homedir(), "Documents", "SPLACE", "iqtree");
    try { fs.rmSync(outputDir, { recursive: true, force: true }); } catch {}
    fs.mkdirSync(outputDir, { recursive: true });

    const nexusFile = path.join(outputDir, "concatenated.nex");
    const partitionFile = path.join(outputDir, "partitions.txt");
    const concatPrefix = path.join(outputDir, "concat_tree");

    fs.writeFileSync(nexusFile, nexus, "utf8");
    fs.writeFileSync(partitionFile, partition, "utf8");

    // Resolve IQ-TREE binary: bundled first, then system PATH in order iqtree3 > iqtree2 > iqtree
    function findIqtree() {
        const names = ["iqtree3", "iqtree2", "iqtree"];
        for (const name of names) {
            const local = resolveBin(name);
            if (local !== name) return local;
        }
        for (const name of names) {
            try {
                const r = spawnSync("which", [name], { encoding: "utf8", shell: false });
                if (r.status === 0 && r.stdout.trim()) return name;
            } catch { }
        }
        return "iqtree3";
    }
    const iqtreeBin = findIqtree();

    // Helper: run one IQ-TREE process and stream output
    function runOne(args, label) {
        return new Promise((resolve) => {
            send("iqtree-progress", { message: `\n[${label}] ${iqtreeBin} ${args.join(" ")}` });
            const proc = spawn(iqtreeBin, args, { shell: true, stdio: ["ignore", "pipe", "pipe"] });
            proc.stdout.on("data", (d) => {
                const msg = d.toString().trim();
                if (msg) send("iqtree-progress", { message: `[${label}] ${msg}` });
            });
            proc.stderr.on("data", (d) => {
                const msg = d.toString().trim();
                if (msg) send("iqtree-progress", { message: `[${label}] ${msg}` });
            });
            proc.on("close", (code) => resolve({ code, label }));
            proc.on("error", (err) => resolve({ code: -1, label, error: err.message }));
        });
    }

    const outgroupArgs = outgroup && outgroup.length ? ["-o", outgroup.join(",")] : [];

    // === Concatenated tree ===
    const concatArgs = [
        "-s", nexusFile,
        "-p", partitionFile,
        "--prefix", concatPrefix,
        "-T", String(threads),
        "--redo",
        ...outgroupArgs,
        ...params,
    ];
    const concatResult = await runOne(concatArgs, "Concatenated");

    // === Per-gene trees (optional) ===
    const geneResults = [];
    if (perGene && geneFiles) {
        const genesDir = path.join(outputDir, "gene_trees");
        fs.mkdirSync(genesDir, { recursive: true });
        const geneNames = Object.keys(geneFiles);
        for (const gene of geneNames) {
            const fastaFile = path.join(genesDir, `${gene}.fasta`);
            const genePrefix = path.join(genesDir, `${gene}_tree`);
            fs.writeFileSync(fastaFile, geneFiles[gene], "utf8");
            const geneArgs = [
                "-s", fastaFile,
                "--prefix", genePrefix,
                "-T", String(threads),
                "--redo",
                ...outgroupArgs,
                ...params,
            ];
            const res = await runOne(geneArgs, gene);
            geneResults.push(res);
        }
    }

    // Collect output files
    let files = [];
    try {
        files = fs.readdirSync(outputDir)
            .filter(f => !f.endsWith(".fasta"))
            .map(f => f);
        if (perGene) {
            const genesDir = path.join(outputDir, "gene_trees");
            try {
                fs.readdirSync(genesDir).forEach(f => files.push(`gene_trees/${f}`));
            } catch { }
        }
    } catch { }

    const success = concatResult.code === 0;
    send("iqtree-done", {
        success,
        outputDir,
        files,
        error: !success ? `Concatenated tree exit code ${concatResult.code}` : null,
        geneResults: geneResults.map(r => ({ gene: r.label, success: r.code === 0 })),
    });
});

// -----------------------------------------------------------------------
// IPC: save ZIP of one or more directories
// { sourceDirs: string[], suggestedName }
// -----------------------------------------------------------------------
ipcMain.handle("save-zip", async (event, { sourceDirs, suggestedName }) => {
    const { filePath, canceled } = await dialog.showSaveDialog({
        defaultPath: path.join(os.homedir(), "Documents", suggestedName || "splace_results.zip"),
        filters: [{ name: "ZIP archive", extensions: ["zip"] }],
    });
    if (canceled || !filePath) return { cancelled: true };

    // Build a staging dir, copy each source into named subfolder, zip it
    const stagingDir = fs.mkdtempSync(path.join(os.tmpdir(), "splace-zip-"));
    try {
        for (const srcDir of sourceDirs) {
            if (!srcDir || !fs.existsSync(srcDir)) continue;
            const folderName = path.basename(srcDir);
            const dest = path.join(stagingDir, folderName);
            fs.mkdirSync(dest, { recursive: true });
            // copy recursively using cp
            const cpResult = spawnSync("cp", ["-r", srcDir + "/.", dest], { shell: false });
            if (cpResult.status !== 0) {
                // fallback: copy files one by one
                for (const f of fs.readdirSync(srcDir, { recursive: true })) {
                    try {
                        const full = path.join(srcDir, f);
                        if (fs.statSync(full).isFile()) {
                            const rel = path.relative(srcDir, full);
                            const target = path.join(dest, rel);
                            fs.mkdirSync(path.dirname(target), { recursive: true });
                            fs.copyFileSync(full, target);
                        }
                    } catch {}
                }
            }
        }
    } catch (e) {
        return { success: false, error: e.message };
    }

    return new Promise((resolve) => {
        const proc = spawn("zip", ["-r", filePath, "."], {
            cwd: stagingDir, shell: true, stdio: ["ignore", "pipe", "pipe"],
        });
        proc.on("close", (code) => {
            try { fs.rmSync(stagingDir, { recursive: true, force: true }); } catch {}
            if (code === 0) resolve({ success: true, filePath });
            else resolve({ success: false, error: `zip exit code ${code}` });
        });
        proc.on("error", (err) => {
            try { fs.rmSync(stagingDir, { recursive: true, force: true }); } catch {}
            resolve({ success: false, error: err.message });
        });
    });
});
