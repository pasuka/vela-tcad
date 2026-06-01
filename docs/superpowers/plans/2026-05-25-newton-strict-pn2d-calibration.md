# Newton Strict pn2d Calibration Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Turn the current pn2d Sentaurus loop from finite diagnostic output into a strict Newton-handoff calibration gate with node-level junction doping handled explicitly.

**Architecture:** Keep `gummel_newton` as a two-stage DC sweep method, but make every accepted pn2d faithful point prove that coupled Newton accepted the final state. Add solver-stage diagnostics to persisted CSV/report artifacts, fix reference comparison so curves are compared at matching biases, add deterministic handling and metadata for compensated TDR junction nodes, then tighten pn2d IV/BV gates. Add reusable reference config coverage for existing checked-in device fixtures after the pn2d gate is honest.

**Tech Stack:** C++20, CMake/Ninja, Catch2, Python stdlib regression tests, existing Sentaurus import tooling, JSON reference configs.

---

## Current Result Review

Observed on 2026-05-25 with:

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
ctest --test-dir build --output-on-failure -R "dc_sweep|sentaurus_sample"
python scripts\sentaurus_import.py reference --config reference_tcad\pn2d\pn2d_reference.json --source-dir reference_tcad\pn2d --output-dir build\pn2d_review --tdr-importer build\sentaurus_import.exe --runner build\vela_example_runner.exe
```

Current status:

- Focused `dc_sweep` and `sentaurus_sample_integration` tests pass.
- Generated pn2d faithful decks preserve `node_doping_file: "doping.csv"` and run finite IV/BV CSVs.
- The pn2d config is not strict Newton: `reference_tcad/pn2d/pn2d_reference.json` sets `handoff.fallback: "gummel_on_newton_failure"`, `require_gummel_convergence: false`, and `newton_max_iter: 0`.
- A strict probe with `fallback: "none"`, `require_gummel_convergence: true`, and no `newton_max_iter: 0` stops at the first IV point and writes only one non-converged row.
- The current CSV does not persist `solver_method`, `gummel_iterations`, `newton_iterations`, or `handoff_stage`, so external reports cannot tell whether a point came from Newton or Gummel fallback.
- The comparison tool currently zips candidate and reference rows by row order, not by bias. For pn2d IV, this compares Vela `0.0..1.0 V` points against the first eleven Sentaurus points clustered near `0.0..0.0326 V`, producing misleading order-of-magnitude numbers.
- In the imported pn2d `doping.csv`, 33 of 882 nodes have both donors and acceptors equal to `1e17 cm^-3`; these are junction/interface compensated nodes produced by merging region-local TDR field values onto a single global node.

These facts explain the next phase: first make solver provenance and curve comparison trustworthy, then remove the Gummel fallback, then calibrate physical trends.

---

## File Structure

- Modify `src/simulation/DCSweep.cpp`
  - Write solver-stage diagnostics into every sweep CSV.
  - Add separate hybrid initializer controls so Gummel and Newton iteration budgets are not tied to one `max_iter`.
- Modify `scripts/compare_reference_curves.py`
  - Align curves by `bias_V`.
  - Support comparison windows and candidate current scale.
- Modify `scripts/sentaurus_import.py`
  - Pass comparison window/scale options from reference config.
  - Add solver-stage summary to manifests.
- Modify `src/io/SentaurusTdrReader.cpp`
  - Detect compensated junction nodes while writing `doping.csv`.
  - Write `doping_metadata.json` with compensation diagnostics.
- Modify `tests/test_sentaurus_tdr_reader.cpp`
  - Cover compensated interface node metadata.
- Modify `tests/test_dc_sweep.cpp`
  - Cover CSV solver-stage diagnostics and hybrid initializer controls.
- Modify `tests/regression/test_reference_tcad_tools.py`
  - Cover bias-aware comparison and reference config discovery.
- Modify `tests/regression/test_sentaurus_import_tools.py`
  - Cover generated strict pn2d config and comparison options.
- Modify `tests/regression/test_sentaurus_sample_integration.py`
  - Require strict Newton handoff for pn2d once solver fixes land.
- Modify `reference_tcad/pn2d/pn2d_reference.json`
  - Remove diagnostic fallback from the faithful pn2d path.
  - Add bias-aware comparison options.
- Create `reference_tcad/nmos2d/nmos2d_reference.json`
- Create `reference_tcad/pmos2d/pmos2d_reference.json`
- Create `reference_tcad/ldmos2d/ldmos2d_reference.json`
- Create `reference_tcad/igbt2d/igbt2d_reference.json`
- Modify `docs/validation/pn2d_sentaurus_comparison.md`
  - Record strict Newton status and calibrated-gate thresholds.
- Modify `docs/config_schema.md`
  - Document persisted solver diagnostics, hybrid initializer controls, and comparison options.

---

## Task 1: Persist Solver Handoff Provenance

**Files:**
- Modify: `src/simulation/DCSweep.cpp`
- Modify: `docs/config_schema.md`
- Test: `tests/test_dc_sweep.cpp`

- [ ] **Step 1: Write the failing CSV diagnostics test**

Add this helper near existing CSV helpers in `tests/test_dc_sweep.cpp`:

```cpp
static std::string readTextFile(const std::filesystem::path& path)
{
    std::ifstream in(path);
    REQUIRE(in.is_open());
    std::ostringstream ss;
    ss << in.rdbuf();
    return ss.str();
}
```

Add this test next to the existing `gummel_newton` tests:

```cpp
TEST_CASE("DCSweep: CSV records hybrid solver handoff provenance",
          "[dc_sweep][gummel_newton][csv]")
{
    const auto dir = makeUniqueSweepDir();
    const ScopedDirectoryCleanup cleanup{dir};
    std::filesystem::create_directories(dir);
    const auto meshPath = writePNMesh(dir);
    const auto csvPath = dir / "handoff_columns.csv";
    const auto cfgPath = writeSweepConfig(dir, meshPath, csvPath, {
        {"start", 0.0},
        {"stop", 0.0},
        {"step", 0.25},
        {"write_vtk", false}
    }, {
        {"method", "gummel_newton"},
        {"max_iter", 12},
        {"reltol", 1.0e-8},
        {"abstol", 1.0e-18},
        {"damping_psi", 0.35},
        {"line_search", true},
        {"warm_start", true},
        {"verbose", false},
        {"handoff", {{"fallback", "none"}}}
    });

    DCSweep sweep;
    const DCSweepResult result = sweep.runWithResult(cfgPath.string());
    REQUIRE(result.points.size() == 1);
    REQUIRE(result.points.front().handoffStage == "newton");

    const std::string csv = readTextFile(csvPath);
    REQUIRE(csv.find("solver_method,gummel_iterations,newton_iterations,handoff_stage") !=
            std::string::npos);
    REQUIRE(csv.find("gummel_newton") != std::string::npos);
    REQUIRE(csv.find(",newton,") != std::string::npos);
}
```

- [ ] **Step 2: Run the focused test and confirm it fails**

Run:

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
cmake --build build --parallel
ctest --test-dir build --output-on-failure -R dc_sweep
```

Expected: the new test fails because the CSV header does not include solver-stage columns.

- [ ] **Step 3: Add solver-stage columns to the CSV header**

In `src/simulation/DCSweep.cpp`, extend the base header from:

```cpp
std::vector<std::string> header = {"mode", "bias_contact", "bias_V",
    "current_contact", "current_electron", "current_hole", "current_total",
    "converged", "iterations", "step_diagnostics", "validation_diagnostics"};
```

to:

```cpp
std::vector<std::string> header = {"mode", "bias_contact", "bias_V",
    "current_contact", "current_electron", "current_hole", "current_total",
    "converged", "iterations", "solver_method", "gummel_iterations",
    "newton_iterations", "handoff_stage", "step_diagnostics",
    "validation_diagnostics"};
```

- [ ] **Step 4: Add solver-stage values to each CSV row**

In `recordPoint`, extend the base `row` construction from:

```cpp
std::vector<std::string> row = {
    toString(sweep.mode),
    sweep.contact,
    formatReal(point.bias),
    sweep.currentContact,
    formatReal(point.electronCurrent),
    formatReal(point.holeCurrent),
    formatReal(point.totalCurrent),
    point.converged ? "1" : "0",
    std::to_string(point.iterations),
    stepDiagnostics(point),
    point.validationDiagnostics};
```

to:

```cpp
std::vector<std::string> row = {
    toString(sweep.mode),
    sweep.contact,
    formatReal(point.bias),
    sweep.currentContact,
    formatReal(point.electronCurrent),
    formatReal(point.holeCurrent),
    formatReal(point.totalCurrent),
    point.converged ? "1" : "0",
    std::to_string(point.iterations),
    point.solverMethod,
    std::to_string(point.gummelIterations),
    std::to_string(point.newtonIterations),
    point.handoffStage,
    stepDiagnostics(point),
    point.validationDiagnostics};
```

- [ ] **Step 5: Document the CSV columns**

In `docs/config_schema.md`, add under DC sweep output:

```markdown
DC sweep CSVs include solver provenance columns:
- `solver_method`: selected nonlinear path, such as `gummel`, `newton`, or `gummel_newton`
- `gummel_iterations`: iterations used by the Gummel stage for this bias point
- `newton_iterations`: iterations used by the coupled Newton stage for this bias point
- `handoff_stage`: final accepted stage or failure stage, such as `newton`, `gummel_failed`, `newton_failed`, or `gummel_fallback`
```

- [ ] **Step 6: Run tests and commit**

Run:

```powershell
cmake --build build --parallel
ctest --test-dir build --output-on-failure -R dc_sweep
```

Expected: all `dc_sweep` tests pass.

Commit:

```powershell
git add src/simulation/DCSweep.cpp docs/config_schema.md tests/test_dc_sweep.cpp
git commit -m "Persist DC sweep solver handoff diagnostics"
```

---

## Task 2: Compare Curves At Matching Bias Points

**Files:**
- Modify: `scripts/compare_reference_curves.py`
- Test: `tests/regression/test_reference_tcad_tools.py`

- [ ] **Step 1: Add failing bias-alignment coverage**

Add this test after `test_compare_reference_curves_enforces_single_curve_thresholds` in `tests/regression/test_reference_tcad_tools.py`:

```python
def test_compare_reference_curves_interpolates_by_bias(self) -> None:
    with tempfile.TemporaryDirectory(prefix="vela_reference_bias_match_") as tmp:
        root = Path(tmp)
        reference = root / "reference.csv"
        candidate = root / "candidate.csv"
        out_json = root / "report.json"
        out_md = root / "report.md"
        self._write_csv(reference, ["bias_V", "current_total"], [
            [0.0, 1.0e-12],
            [0.5, 1.0e-9],
            [1.0, 1.0e-6],
        ])
        self._write_csv(candidate, ["bias_V", "current_total"], [
            [0.0, -1.0e-12],
            [0.25, -3.0e-11],
            [0.75, -3.0e-8],
            [1.0, -1.0e-6],
        ])

        subprocess.run([
            sys.executable,
            str(REPO / "scripts" / "compare_reference_curves.py"),
            "--reference", str(reference),
            "--candidate", str(candidate),
            "--output-json", str(out_json),
            "--output-md", str(out_md),
            "--kind", "iv",
            "--candidate-scale", "-1.0",
            "--bias-min", "0.5",
            "--bias-max", "1.0",
            "--max-orders-of-magnitude", "0.25",
            "--require-trend-match",
        ], check=True, cwd=REPO)

        report = json.loads(out_json.read_text())
        self.assertEqual(report["status"], "pass")
        self.assertEqual(report["iv"]["points_compared"], 2)
        self.assertEqual(report["iv"]["reference_bias_range"], [0.5, 1.0])
        self.assertEqual(report["iv"]["candidate_scale"], -1.0)
```

- [ ] **Step 2: Run and confirm failure**

Run:

```powershell
python -m unittest tests.regression.test_reference_tcad_tools -v
```

Expected: failure because the compare script has no `--candidate-scale`, `--bias-min`, or `--bias-max` arguments and compares rows by position.

- [ ] **Step 3: Add interpolation helpers**

In `scripts/compare_reference_curves.py`, add:

```python
def finite_pairs(rows: list[dict[str, str]], value_column: str, scale: float) -> list[tuple[float, float]]:
    pairs: list[tuple[float, float]] = []
    for row in rows:
        if "bias_V" not in row:
            continue
        bias = float(row["bias_V"])
        value = float(row[value_column]) * scale
        if math.isfinite(bias) and math.isfinite(value):
            pairs.append((bias, value))
    return sorted(pairs)


def interpolate_at(pairs: list[tuple[float, float]], bias: float) -> float | None:
    if not pairs:
        return None
    if bias < pairs[0][0] or bias > pairs[-1][0]:
        return None
    for existing_bias, value in pairs:
        if abs(existing_bias - bias) <= max(abs(bias), 1.0) * 1.0e-12:
            return value
    for (b0, v0), (b1, v1) in zip(pairs, pairs[1:]):
        if b0 <= bias <= b1 and b1 != b0:
            t = (bias - b0) / (b1 - b0)
            return v0 + t * (v1 - v0)
    return None


def aligned_values(reference_rows: list[dict[str, str]],
                   candidate_rows: list[dict[str, str]],
                   ref_col: str,
                   cand_col: str,
                   candidate_scale: float,
                   bias_min: float | None,
                   bias_max: float | None) -> tuple[list[float], list[float], list[float]]:
    ref_pairs = finite_pairs(reference_rows, ref_col, 1.0)
    cand_pairs = finite_pairs(candidate_rows, cand_col, candidate_scale)
    biases: list[float] = []
    ref_values: list[float] = []
    cand_values: list[float] = []
    for bias, ref_value in ref_pairs:
        if bias_min is not None and bias < bias_min:
            continue
        if bias_max is not None and bias > bias_max:
            continue
        cand_value = interpolate_at(cand_pairs, bias)
        if cand_value is None:
            continue
        biases.append(bias)
        ref_values.append(ref_value)
        cand_values.append(cand_value)
    return biases, ref_values, cand_values
```

- [ ] **Step 4: Wire comparison options**

Change `compare_series` to accept `candidate_scale`, `bias_min`, and `bias_max`, and replace the old `values(...)` calls with:

```python
biases, ref_values, cand_values = aligned_values(
    reference_rows,
    candidate_rows,
    ref_col,
    cand_col,
    candidate_scale,
    bias_min,
    bias_max,
)
```

Add these fields to the returned dict:

```python
"reference_bias_range": [biases[0], biases[-1]] if biases else None,
"candidate_scale": candidate_scale,
```

Add parser options:

```python
parser.add_argument("--candidate-scale", type=float, default=1.0)
parser.add_argument("--bias-min", type=float)
parser.add_argument("--bias-max", type=float)
```

Pass those arguments into each `compare_series(...)` call.

- [ ] **Step 5: Run tests and commit**

Run:

```powershell
python -m unittest tests.regression.test_reference_tcad_tools -v
```

Expected: all reference tool tests pass.

Commit:

```powershell
git add scripts/compare_reference_curves.py tests/regression/test_reference_tcad_tools.py
git commit -m "Compare reference curves by bias"
```

---

## Task 3: Report And Resolve Compensated TDR Junction Nodes

**Files:**
- Modify: `src/io/SentaurusTdrReader.cpp`
- Test: `tests/test_sentaurus_tdr_reader.cpp`

- [ ] **Step 1: Add failing metadata coverage**

Extend the synthetic TDR export test in `tests/test_sentaurus_tdr_reader.cpp` so the fixture has one shared node with equal donor and acceptor active concentration. After export, assert:

```cpp
const auto metadata = nlohmann::json::parse(readFile(outDir / "doping_metadata.json"));
REQUIRE(metadata["compensated_nodes"]["count"].get<int>() == 1);
REQUIRE(metadata["compensated_nodes"]["nodes"][0]["node_id"].get<int>() == 1);
REQUIRE(metadata["compensated_nodes"]["nodes"][0]["donors_cm3"].get<double>() == Catch::Approx(1.0e17));
REQUIRE(metadata["compensated_nodes"]["nodes"][0]["acceptors_cm3"].get<double>() == Catch::Approx(1.0e17));
REQUIRE(metadata["compensated_nodes"]["nodes"][0]["policy"].get<std::string>() == "reported");
```

- [ ] **Step 2: Run and confirm failure**

Run:

```powershell
cmake --build build --parallel
ctest --test-dir build --output-on-failure -R sentaurus_tdr_reader
```

Expected: failure because `doping_metadata.json` is not written.

- [ ] **Step 3: Write compensation metadata**

In `src/io/SentaurusTdrReader.cpp`, after writing `doping.csv`, add:

```cpp
nlohmann::json dopingMetadata;
dopingMetadata["schema"] = "vela.sentaurus_tdr.doping_metadata.v1";
dopingMetadata["compensated_nodes"] = {
    {"policy", "reported"},
    {"count", 0},
    {"nodes", nlohmann::json::array()},
};
for (std::size_t i = 0; i < inventory.vertices.size(); ++i) {
    const double donor = donors[i];
    const double acceptor = acceptors[i];
    const double scale = std::max({std::abs(donor), std::abs(acceptor), 1.0});
    if (donor > 0.0 && acceptor > 0.0 &&
        std::abs(donor - acceptor) <= 1.0e-6 * scale) {
        dopingMetadata["compensated_nodes"]["nodes"].push_back({
            {"node_id", i},
            {"donors_cm3", donor},
            {"acceptors_cm3", acceptor},
            {"policy", "reported"},
        });
    }
}
dopingMetadata["compensated_nodes"]["count"] =
    dopingMetadata["compensated_nodes"]["nodes"].size();
{
    std::ofstream out(outDir / "doping_metadata.json");
    out << dopingMetadata.dump(2) << "\n";
}
```

- [ ] **Step 4: Run tests and commit**

Run:

```powershell
cmake --build build --parallel
ctest --test-dir build --output-on-failure -R sentaurus_tdr_reader
```

Expected: the Sentaurus TDR reader tests pass and the pn2d import reports 33 compensated nodes in `doping_metadata.json`.

Commit:

```powershell
git add src/io/SentaurusTdrReader.cpp tests/test_sentaurus_tdr_reader.cpp
git commit -m "Report compensated Sentaurus junction doping nodes"
```

---

## Task 4: Add A Deterministic Junction Compensation Policy

**Files:**
- Modify: `src/io/SentaurusTdrReader.cpp`
- Modify: `scripts/sentaurus_import.py`
- Modify: `docs/config_schema.md`
- Test: `tests/test_sentaurus_tdr_reader.cpp`
- Test: `tests/regression/test_sentaurus_import_tools.py`

- [ ] **Step 1: Add failing policy tests**

In `tests/test_sentaurus_tdr_reader.cpp`, add a synthetic TDR export case with a shared junction node and signed `DopingConcentration` fields on two material regions. Run the importer with:

```cpp
SentaurusTdrExportOptions options;
options.compensatedDopingPolicy = "dominant_signed_region";
exportSentaurusTdrToReferenceCsv(tdrPath, outDir, options);
```

Assert the exported shared node is no longer exactly compensated:

```cpp
const std::string doping = readFile(outDir / "doping.csv");
REQUIRE(doping.find("1,1e+17,0") != std::string::npos);
const auto metadata = nlohmann::json::parse(readFile(outDir / "doping_metadata.json"));
REQUIRE(metadata["compensated_nodes"]["nodes"][0]["policy"].get<std::string>() ==
        "dominant_signed_region");
REQUIRE(metadata["compensated_nodes"]["nodes"][0]["resolved"].get<bool>());
```

In `tests/regression/test_sentaurus_import_tools.py`, add a config assertion:

```python
self.assertEqual(config["tdr_doping"]["compensated_node_policy"], "dominant_signed_region")
```

and assert the generated manifest includes:

```python
self.assertIn("doping_metadata.json", summary["generated"])
```

- [ ] **Step 2: Run and confirm failure**

Run:

```powershell
ctest --test-dir build --output-on-failure -R sentaurus_tdr_reader
python -m unittest tests.regression.test_sentaurus_import_tools -v
```

Expected: failures because export options and config plumbing do not exist.

- [ ] **Step 3: Add export options**

In the Sentaurus TDR reader public header, add:

```cpp
struct SentaurusTdrExportOptions {
    std::string compensatedDopingPolicy = "reported";
};
```

Add an overload:

```cpp
void exportSentaurusTdrToReferenceCsv(const std::filesystem::path& filename,
                                      const std::filesystem::path& outDir,
                                      const SentaurusTdrExportOptions& options);
```

Keep the existing two-argument function and implement it as:

```cpp
exportSentaurusTdrToReferenceCsv(filename, outDir, SentaurusTdrExportOptions{});
```

- [ ] **Step 4: Implement dominant signed resolution**

While collecting TDR fields, keep signed aggregate `DopingConcentration` contributions per node:

```cpp
std::vector<std::vector<double>> signedDopingByNode(inventory.vertices.size());
```

When `field.name == "DopingConcentration"`, push each region-local value into the node vector. Before writing `doping.csv`, resolve compensated nodes when requested:

```cpp
if (options.compensatedDopingPolicy == "dominant_signed_region") {
    for (std::size_t i = 0; i < inventory.vertices.size(); ++i) {
        const double donor = donors[i];
        const double acceptor = acceptors[i];
        const double scale = std::max({std::abs(donor), std::abs(acceptor), 1.0});
        if (!(donor > 0.0 && acceptor > 0.0 &&
              std::abs(donor - acceptor) <= 1.0e-6 * scale)) {
            continue;
        }
        double signedPick = 0.0;
        for (double value : signedDopingByNode[i]) {
            if (std::abs(value) > std::abs(signedPick))
                signedPick = value;
        }
        if (signedPick > 0.0) {
            donors[i] = std::abs(signedPick);
            acceptors[i] = 0.0;
        } else if (signedPick < 0.0) {
            donors[i] = 0.0;
            acceptors[i] = std::abs(signedPick);
        }
    }
} else if (options.compensatedDopingPolicy != "reported") {
    throw std::invalid_argument(
        "SentaurusTdrExportOptions: compensatedDopingPolicy must be 'reported' or "
        "'dominant_signed_region'.");
}
```

Record both original and resolved values in `doping_metadata.json`.

- [ ] **Step 5: Wire the Python reference config**

In `scripts/sentaurus_import.py`, read:

```python
tdr_doping = config.get("tdr_doping", {})
compensated_policy = tdr_doping.get("compensated_node_policy", "reported")
```

Pass the policy to the C++ importer command as:

```python
cmd.extend(["--compensated-doping-policy", compensated_policy])
```

Add the same CLI argument to `src/tools/sentaurus_import.cpp` and populate `SentaurusTdrExportOptions`.

- [ ] **Step 6: Document config**

In `docs/config_schema.md`, add this text:

````markdown
Sentaurus reference import may include:

```json
"tdr_doping": {
  "compensated_node_policy": "reported"
}
```

Supported policies:
- `reported`: preserve `doping.csv` exactly as merged from region-local TDR fields and report compensated nodes in `doping_metadata.json`.
- `dominant_signed_region`: when a global node receives equal donor and acceptor active concentrations, use the signed `DopingConcentration` field with the largest magnitude to choose a single majority dopant for the node, and record the rewrite in `doping_metadata.json`.
````

- [ ] **Step 7: Run tests and commit**

Run:

```powershell
cmake --build build --parallel
ctest --test-dir build --output-on-failure -R sentaurus_tdr_reader
python -m unittest tests.regression.test_sentaurus_import_tools -v
```

Expected: tests pass and pn2d import can opt into `dominant_signed_region`.

Commit:

```powershell
git add src/io/SentaurusTdrReader.cpp include/vela/io/SentaurusTdrReader.h src/tools/sentaurus_import.cpp scripts/sentaurus_import.py docs/config_schema.md tests/test_sentaurus_tdr_reader.cpp tests/regression/test_sentaurus_import_tools.py
git commit -m "Resolve compensated TDR junction doping nodes"
```

---

## Task 5: Separate Gummel Initializer And Newton Budgets

**Files:**
- Modify: `src/simulation/DCSweep.cpp`
- Modify: `docs/config_schema.md`
- Test: `tests/test_dc_sweep.cpp`

- [ ] **Step 1: Add failing budget test**

Add this test near the hybrid fallback tests:

```cpp
TEST_CASE("DCSweep: hybrid handoff has separate Gummel and Newton iteration budgets",
          "[dc_sweep][gummel_newton]")
{
    const auto dir = makeUniqueSweepDir();
    const ScopedDirectoryCleanup cleanup{dir};
    std::filesystem::create_directories(dir);
    const auto meshPath = writePNMesh(dir);
    const auto csvPath = dir / "hybrid_budget.csv";
    const auto cfgPath = writeSweepConfig(dir, meshPath, csvPath, {
        {"start", 0.0},
        {"stop", 0.0},
        {"step", 0.25},
        {"write_vtk", false}
    }, {
        {"method", "gummel_newton"},
        {"max_iter", 1},
        {"reltol", 1.0e-8},
        {"abstol", 1.0e-18},
        {"damping_psi", 0.35},
        {"line_search", true},
        {"warm_start", true},
        {"verbose", false},
        {"handoff", {
            {"fallback", "none"},
            {"gummel_max_iter", 20},
            {"newton_max_iter", 12}
        }}
    });

    DCSweep sweep;
    const DCSweepResult result = sweep.runWithResult(cfgPath.string());
    REQUIRE(result.points.size() == 1);
    REQUIRE(result.points.front().converged);
    REQUIRE(result.points.front().gummelIterations > 1);
    REQUIRE(result.points.front().newtonIterations <= 12);
    REQUIRE(result.points.front().handoffStage == "newton");
}
```

- [ ] **Step 2: Run and confirm failure**

Run:

```powershell
ctest --test-dir build --output-on-failure -R dc_sweep
```

Expected: failure because `solver.handoff.gummel_max_iter` is ignored.

- [ ] **Step 3: Add `gummel_max_iter` parsing**

Extend `HybridHandoffConfig`:

```cpp
int gummelMaxIter = -1;
```

In `hybridHandoffConfigFromJson`, add:

```cpp
if (handoff.contains("gummel_max_iter")) {
    hybrid.gummelMaxIter = handoff.at("gummel_max_iter").get<int>();
    if (hybrid.gummelMaxIter < 0)
        throw std::invalid_argument(
            "DCSweep: solver.handoff.gummel_max_iter must be non-negative.");
}
```

- [ ] **Step 4: Apply the initializer budget**

In the `SolverMethod::GummelNewton` branch before `runGummel(...)`, add:

```cpp
GummelConfig initializerGummel = gummel;
if (hybrid.gummelMaxIter >= 0)
    initializerGummel.maxIter = hybrid.gummelMaxIter;
```

Use `initializerGummel` in both `runGummel(...)` calls for the hybrid branch.

- [ ] **Step 5: Document config**

In `docs/config_schema.md`, extend hybrid keys:

```markdown
- `handoff.gummel_max_iter`: optional non-negative integer overriding only the Gummel initializer iteration limit.
- `handoff.newton_max_iter`: optional non-negative integer overriding only the Newton handoff iteration limit.
```

- [ ] **Step 6: Run tests and commit**

Run:

```powershell
cmake --build build --parallel
ctest --test-dir build --output-on-failure -R dc_sweep
```

Expected: all `dc_sweep` tests pass.

Commit:

```powershell
git add src/simulation/DCSweep.cpp docs/config_schema.md tests/test_dc_sweep.cpp
git commit -m "Separate hybrid initializer and Newton budgets"
```

---

## Task 6: Make pn2d Faithful Deck Strict Newton

**Files:**
- Modify: `reference_tcad/pn2d/pn2d_reference.json`
- Modify: `scripts/sentaurus_import.py`
- Modify: `docs/validation/pn2d_sentaurus_comparison.md`
- Test: `tests/regression/test_sentaurus_import_tools.py`
- Test: `tests/regression/test_sentaurus_sample_integration.py`

- [ ] **Step 1: Add failing strict assertions**

In `tests/regression/test_sentaurus_sample_integration.py`, change the pn2d deck assertions to:

```python
self.assertEqual(iv_deck["solver"]["handoff"]["fallback"], "none")
self.assertTrue(iv_deck["solver"]["handoff"]["require_gummel_convergence"])
self.assertGreater(iv_deck["solver"]["handoff"]["gummel_max_iter"], 40)
self.assertGreater(iv_deck["solver"]["handoff"]["newton_max_iter"], 0)
```

After reading `faithful_iv`, add:

```python
for row in self._read_curve(faithful_iv):
    self.assertEqual(row["solver_method"], "gummel_newton")
    self.assertEqual(row["handoff_stage"], "newton")
    self.assertGreater(int(row["newton_iterations"]), 0)
```

Add the same `handoff_stage == "newton"` assertion for `faithful_bv`.

- [ ] **Step 2: Run and confirm failure**

Run:

```powershell
python -m unittest tests.regression.test_sentaurus_sample_integration -v
```

Expected: failure because pn2d currently uses diagnostic fallback and CSVs do not yet show all rows ending in `newton`.

- [ ] **Step 3: Tighten pn2d config**

Change `reference_tcad/pn2d/pn2d_reference.json` top-level `vela_solver` to:

```json
"vela_solver": {
  "method": "gummel_newton",
  "max_iter": 30,
  "reltol": 1.0e-8,
  "abstol": 1.0e-18,
  "damping_psi": 0.2,
  "damping_factor": 1.0,
  "line_search": true,
  "warm_start": true,
  "verbose": false,
  "handoff": {
    "fallback": "none",
    "require_gummel_convergence": true,
    "gummel_max_iter": 120,
    "newton_max_iter": 30
  }
},
"tdr_doping": {
  "compensated_node_policy": "dominant_signed_region"
}
```

- [ ] **Step 4: Add comparison options to pn2d simulations**

For IV, add:

```json
"comparison": {
  "candidate_scale": -1.0,
  "bias_min": 0.2,
  "bias_max": 1.0,
  "max_orders_of_magnitude": 2.0,
  "require_trend_match": true,
  "min_points": 6
}
```

For BV, add:

```json
"comparison": {
  "candidate_scale": 1.0,
  "bias_min": 5.0,
  "bias_max": 50.0,
  "max_orders_of_magnitude": 6.0,
  "require_trend_match": false,
  "min_points": 8
}
```

- [ ] **Step 5: Pass comparison options to the compare script**

In `scripts/sentaurus_import.py`, when building `compare_cmd`, read:

```python
comparison = sim.get("comparison", {})
```

Append options when present:

```python
if "candidate_scale" in comparison:
    compare_cmd.extend(["--candidate-scale", str(comparison["candidate_scale"])])
if "bias_min" in comparison:
    compare_cmd.extend(["--bias-min", str(comparison["bias_min"])])
if "bias_max" in comparison:
    compare_cmd.extend(["--bias-max", str(comparison["bias_max"])])
if "max_orders_of_magnitude" in comparison:
    compare_cmd.extend([
        "--max-orders-of-magnitude",
        str(comparison["max_orders_of_magnitude"]),
    ])
if comparison.get("require_trend_match", sim.get("require_trend_match", False)):
    compare_cmd.append("--require-trend-match")
if "min_points" in comparison:
    compare_cmd.extend(["--min-points", str(comparison["min_points"])])
```

- [ ] **Step 6: Run the pn2d import gate**

Run:

```powershell
python scripts\sentaurus_import.py reference --config reference_tcad\pn2d\pn2d_reference.json --source-dir reference_tcad\pn2d --output-dir build\pn2d_strict_gate --tdr-importer build\sentaurus_import.exe --runner build\vela_example_runner.exe
```

Expected after Tasks 1-5: `build\pn2d_strict_gate\vela\pn2d_iv.csv` and `pn2d_bv.csv` exist, every converged row has `handoff_stage,newton`, and comparison reports are generated with bias-aligned metrics.

- [ ] **Step 7: Update validation doc**

In `docs/validation/pn2d_sentaurus_comparison.md`, replace the diagnostic fallback paragraph with:

```markdown
Faithful pn2d decks use `solver.method: "gummel_newton"` with strict handoff:
Gummel initializes each point, but accepted faithful IV/BV rows must end with
`handoff_stage: "newton"` and `newton_iterations > 0`. Gummel fallback is no
longer part of the default pn2d gate.

The comparison report aligns curves by `bias_V`; IV uses `candidate_scale:
-1.0` to match the Sentaurus terminal-current orientation and checks the
0.2-1.0 V forward-bias window. BV remains a diagnostic high-field/current
comparison over 5-50 V.
```

- [ ] **Step 8: Run tests and commit**

Run:

```powershell
cmake --build build --parallel
ctest --test-dir build --output-on-failure -R "dc_sweep|sentaurus_sample|sentaurus_tdr_reader"
python -m unittest tests.regression.test_sentaurus_import_tools -v
```

Expected: all selected tests pass.

Commit:

```powershell
git add reference_tcad/pn2d/pn2d_reference.json scripts/sentaurus_import.py docs/validation/pn2d_sentaurus_comparison.md tests/regression/test_sentaurus_import_tools.py tests/regression/test_sentaurus_sample_integration.py
git commit -m "Require strict Newton handoff for pn2d"
```

---

## Task 7: Add Reference Configs For Checked-In Device Fixtures

**Files:**
- Create: `reference_tcad/nmos2d/nmos2d_reference.json`
- Create: `reference_tcad/pmos2d/pmos2d_reference.json`
- Create: `reference_tcad/ldmos2d/ldmos2d_reference.json`
- Create: `reference_tcad/igbt2d/igbt2d_reference.json`
- Modify: `tests/regression/test_reference_tcad_tools.py`
- Modify: `reference_tcad/README.md`

- [ ] **Step 1: Add failing config inventory test**

In `tests/regression/test_reference_tcad_tools.py`, add:

```python
def test_checked_in_reference_configs_cover_device_fixtures(self) -> None:
    expected = {
        "nmos2d": ["idvd", "idvg", "idvg_surface", "cv", "bv"],
        "pmos2d": ["idvd", "idvg", "idvg_surface", "cv", "bv"],
        "ldmos2d": ["iv", "bv", "fieldplate"],
        "igbt2d": ["iv", "high_injection_iv", "charge_cv", "bv", "bv_ii"],
    }
    for device, simulations in expected.items():
        with self.subTest(device=device):
            path = REPO / "reference_tcad" / device / f"{device}_reference.json"
            self.assertTrue(path.is_file(), path)
            config = json.loads(path.read_text())
            self.assertEqual(config["case"], device)
            self.assertEqual(config["schema"], "vela.reference_tcad.checked_in.v1")
            self.assertEqual(
                [sim["name"] for sim in config["simulations"]],
                simulations,
            )
            for sim in config["simulations"]:
                vela = REPO / "reference_tcad" / device / "vela" / sim["deck"]
                reference = REPO / "reference_tcad" / device / "reference_curves" / sim["reference_curve"]
                report = REPO / "reference_tcad" / device / "reports" / sim["report_json"]
                self.assertTrue(vela.is_file(), vela)
                self.assertTrue(reference.is_file(), reference)
                self.assertTrue(report.is_file(), report)
```

- [ ] **Step 2: Run and confirm failure**

Run:

```powershell
python -m unittest tests.regression.test_reference_tcad_tools -v
```

Expected: failure because the new reference config files do not exist.

- [ ] **Step 3: Create NMOS config**

Create `reference_tcad/nmos2d/nmos2d_reference.json`:

```json
{
  "schema": "vela.reference_tcad.checked_in.v1",
  "case": "nmos2d",
  "device": "nmos2d",
  "mesh": "vela/mesh.json",
  "simulations": [
    {"name": "idvd", "deck": "simulation_idvd.json", "candidate": "nmos2d_idvd.csv", "reference_curve": "nmos2d_idvd_reference.csv", "report_json": "nmos2d_idvd_comparison.json", "kind": "iv"},
    {"name": "idvg", "deck": "simulation_idvg.json", "candidate": "nmos2d_idvg.csv", "reference_curve": "nmos2d_idvg_reference.csv", "report_json": "nmos2d_idvg_comparison.json", "kind": "iv"},
    {"name": "idvg_surface", "deck": "simulation_idvg_surface.json", "candidate": "nmos2d_idvg_surface.csv", "reference_curve": "nmos2d_idvg_reference.csv", "report_json": "nmos2d_idvg_comparison.json", "kind": "iv"},
    {"name": "cv", "deck": "simulation_cv.json", "candidate": "nmos2d_cv.csv", "reference_curve": "nmos2d_cv_reference.csv", "report_json": "nmos2d_cv_comparison.json", "kind": "cv"},
    {"name": "bv", "deck": "simulation_bv.json", "candidate": "nmos2d_bv.csv", "reference_curve": "nmos2d_bv_reference.csv", "report_json": "nmos2d_bv_comparison.json", "kind": "bv"}
  ]
}
```

- [ ] **Step 4: Create PMOS, LDMOS, and IGBT configs**

Create equivalent files using the existing candidate and report names from each directory:

```json
{
  "schema": "vela.reference_tcad.checked_in.v1",
  "case": "pmos2d",
  "device": "pmos2d",
  "mesh": "vela/mesh.json",
  "simulations": [
    {"name": "idvd", "deck": "simulation_idvd.json", "candidate": "pmos2d_idvd.csv", "reference_curve": "pmos2d_idvd_reference.csv", "report_json": "pmos2d_idvd_comparison.json", "kind": "iv"},
    {"name": "idvg", "deck": "simulation_idvg.json", "candidate": "pmos2d_idvg.csv", "reference_curve": "pmos2d_idvg_reference.csv", "report_json": "pmos2d_idvg_comparison.json", "kind": "iv"},
    {"name": "idvg_surface", "deck": "simulation_idvg_surface.json", "candidate": "pmos2d_idvg_surface.csv", "reference_curve": "pmos2d_idvg_reference.csv", "report_json": "pmos2d_idvg_comparison.json", "kind": "iv"},
    {"name": "cv", "deck": "simulation_cv.json", "candidate": "pmos2d_cv.csv", "reference_curve": "pmos2d_cv_reference.csv", "report_json": "pmos2d_cv_comparison.json", "kind": "cv"},
    {"name": "bv", "deck": "simulation_bv.json", "candidate": "pmos2d_bv.csv", "reference_curve": "pmos2d_bv_reference.csv", "report_json": "pmos2d_bv_comparison.json", "kind": "bv"}
  ]
}
```

```json
{
  "schema": "vela.reference_tcad.checked_in.v1",
  "case": "ldmos2d",
  "device": "ldmos2d",
  "mesh": "vela/mesh.json",
  "simulations": [
    {"name": "iv", "deck": "simulation_iv.json", "candidate": "ldmos2d_iv.csv", "reference_curve": "ldmos2d_iv_reference.csv", "report_json": "ldmos2d_iv_comparison.json", "kind": "iv"},
    {"name": "bv", "deck": "simulation_bv.json", "candidate": "ldmos2d_bv.csv", "reference_curve": "ldmos2d_bv_reference.csv", "report_json": "ldmos2d_bv_comparison.json", "kind": "bv"},
    {"name": "fieldplate", "deck": "simulation_bv_fieldplate.json", "candidate": "ldmos2d_fieldplate.csv", "reference_curve": "ldmos2d_fieldplate_reference.csv", "report_json": "ldmos2d_fieldplate_comparison.json", "kind": "bv"}
  ]
}
```

```json
{
  "schema": "vela.reference_tcad.checked_in.v1",
  "case": "igbt2d",
  "device": "igbt2d",
  "mesh": "vela/mesh.json",
  "simulations": [
    {"name": "iv", "deck": "simulation_iv.json", "candidate": "igbt2d_iv.csv", "reference_curve": "igbt2d_iv_reference.csv", "report_json": "igbt2d_iv_comparison.json", "kind": "iv"},
    {"name": "high_injection_iv", "deck": "simulation_high_injection_iv.json", "candidate": "igbt2d_high_injection_iv.csv", "reference_curve": "igbt2d_high_injection_iv_reference.csv", "report_json": "igbt2d_high_injection_iv_comparison.json", "kind": "iv"},
    {"name": "charge_cv", "deck": "simulation_charge_cv.json", "candidate": "igbt2d_charge_cv.csv", "reference_curve": "igbt2d_charge_cv_reference.csv", "report_json": "igbt2d_charge_cv_comparison.json", "kind": "cv"},
    {"name": "bv", "deck": "simulation_bv.json", "candidate": "igbt2d_bv.csv", "reference_curve": "igbt2d_bv_reference.csv", "report_json": "igbt2d_bv_comparison.json", "kind": "bv"},
    {"name": "bv_ii", "deck": "simulation_bv_ii.json", "candidate": "igbt2d_bv_ii.csv", "reference_curve": "igbt2d_bv_ii_reference.csv", "report_json": "igbt2d_bv_ii_comparison.json", "kind": "bv"}
  ]
}
```

- [ ] **Step 5: Document reference configs**

In `reference_tcad/README.md`, add:

```markdown
Each checked-in fixture may include `<device>_reference.json` with schema
`vela.reference_tcad.checked_in.v1`. These configs inventory the mesh, Vela
decks, candidate CSVs, reference curves, comparison reports, and curve kind for
each reusable sample. They are intentionally metadata-only for checked-in CSV
fixtures; generated Sentaurus imports use `vela.reference_tcad.sentaurus_reference.v1`.
```

- [ ] **Step 6: Run tests and commit**

Run:

```powershell
python -m unittest tests.regression.test_reference_tcad_tools -v
```

Expected: reference TCAD regression tests pass and all four new config files are covered.

Commit:

```powershell
git add reference_tcad/nmos2d/nmos2d_reference.json reference_tcad/pmos2d/pmos2d_reference.json reference_tcad/ldmos2d/ldmos2d_reference.json reference_tcad/igbt2d/igbt2d_reference.json reference_tcad/README.md tests/regression/test_reference_tcad_tools.py
git commit -m "Add reference configs for checked-in device fixtures"
```

---

## Final Verification

Run:

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
cmake --build build --parallel
ctest --test-dir build --output-on-failure
python scripts\sentaurus_import.py reference --config reference_tcad\pn2d\pn2d_reference.json --source-dir reference_tcad\pn2d --output-dir build\pn2d_strict_final --tdr-importer build\sentaurus_import.exe --runner build\vela_example_runner.exe
```

Expected final state:

- Full CTest passes.
- pn2d faithful IV/BV CSVs include solver provenance columns.
- Every accepted pn2d faithful row has `handoff_stage` equal to `newton`.
- pn2d comparison reports align by `bias_V`.
- pn2d IV comparison uses the configured current sign scale and forward-bias comparison window.
- `doping_metadata.json` records compensated junction nodes and the applied policy.
- `reference_tcad/nmos2d`, `pmos2d`, `ldmos2d`, and `igbt2d` each have a checked-in reference config.

---

## Self-Review

- Spec coverage: The plan addresses strict Newton takeover, pn2d IV/BV magnitude/trend comparison, TDR node-level compensated junction doping, and more device sample reference configs.
- Placeholder scan: No step uses placeholder wording or an undefined future component. Each code-changing task has exact files, test commands, expected failure, and implementation snippets.
- Type consistency: `solver_method`, `gummel_iterations`, `newton_iterations`, `handoff_stage`, `gummel_max_iter`, `newton_max_iter`, `candidate_scale`, `bias_min`, `bias_max`, and `compensated_node_policy` are used consistently across tests, scripts, JSON, docs, and C++.
