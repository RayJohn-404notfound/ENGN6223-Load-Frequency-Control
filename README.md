# Single-Area Load Frequency Control

An ENGN6223 control-systems project comparing PI and filtered PID (PIDF) control on a continuous-time power-system benchmark.

**Research question:** Why can controllers with the same crossover frequency and phase margin recover differently after a load increase?

The study connects time responses, loop sensitivity, Nyquist stability, modal analysis, parameter variation, anti-windup and measurement-noise gain. It uses a literature model; it does not establish physical-grid performance or general PIDF superiority.

## Main result

The original PI/PIDF comparison uses a crossover of **2.4 rad/s** and a **65°** phase target. Nominal recovery is **4.655 s / 2.991 s**, but three parameter cases reverse that ordering.

A separate design study selects a **75° PIDF candidate** for further development. It passes **42 of 43** deterministic cases; case 31 misses the 45° phase-margin criterion by **0.842°**. Its noise-to-command gain at 100 rad/s is **70.9% higher** than the original PIDF. No evaluated design passes every original criterion.

See [results and interpretation](docs/RESULTS.md) for the comparison, assumptions and remaining limitations.

## Requirements

- MATLAB and Control System Toolbox. The original numerical study used MATLAB R2026a Update 3.
- Simulink is optional for the explicit-block cross-check and model diagrams.
- Python 3, NumPy and Matplotlib for the final figures. The two preceding Python analyses use the standard library.

## Run the study

Use a working copy: execution replaces generated outputs under `results/` and may regenerate models under `model/`.

**1. In MATLAB**, open the repository folder and run:

```matlab
run_project
```

This runs the original comparison, numerical checks, recovery diagnostics and design grid. The Simulink cross-check runs when Simulink is available. Keep `results/design_sensitivity/protocol.json`: it is an input, despite being under `results/`.

**2. In a terminal at the repository root**, run the saved-result analyses in this order:

```sh
python3 -m pip install -r requirements.txt
python3 results/design_decision_review/review_saved_designs.py
python3 results/report_evidence_summary/audit_saved_evidence.py
python3 results/report_compact/make_report_figures.py
```

The Python scripts reclassify saved metrics, check summary calculations and create four figure pairs. They do not run new dynamical simulations or retune controllers.

## Reproduction check

On 30 October 2026, the full MATLAB entry point completed in a separate copy containing only source files and the design protocol. MATLAB R2026a Update 3 regenerated the results, including 26 feasible designs, 1,118 design-case evaluations and 15 matched pairs. The nominal transfer-function, state-space, ODE and Simulink cross-checks passed. All three Python post-processing scripts then completed using Python 3.10.9, NumPy 1.23.5 and Matplotlib 3.7.0. The development selection remained `PIDF_w24_pm75_N10`, with 42/43 cases meeting all declared criteria.

The committed CSV/JSON snapshots retain the original reference results. This local reproduction checks numerical consistency, not performance on a physical power system.

## Repository guide

| Location | Purpose |
|---|---|
| `run_project.m` | Main MATLAB entry point |
| `matlab/` | Model, controller design, simulations and numerical checks |
| `model/` | Generated explicit-block Simulink model; see its README |
| `results/` | Selected numerical snapshots, the retained protocol and post-processing scripts |
| [`docs/RESULTS.md`](docs/RESULTS.md) | Key findings, evidence map and interpretation limits |

The main MATLAB script also contains supplementary pole-placement, LQR and observer demonstrations. These illustrate course concepts and are outside the main PI/PIDF comparison. The observer example assumes a known zero disturbance.

The assessment report, slides, video, teaching materials and working drafts are kept separately. Saved results are reference outputs; rerunning the code and physical-system validation are distinct checks.

## Literature

- Anwar and Pan (2015), *A new PID load frequency controller design method in frequency domain through direct synthesis approach*. [DOI: 10.1016/j.ijepes.2014.12.024](https://doi.org/10.1016/j.ijepes.2014.12.024). Example 1 supplies the benchmark parameters; its tuning method is not reproduced here.
- Hote and Jain (2018), *PID controller design for load frequency control: Past, present and future challenges*. [DOI: 10.1016/j.ifacol.2018.06.162](https://doi.org/10.1016/j.ifacol.2018.06.162). Background on LFC constraints and validation challenges.

## Project team

Rui Jiang, Zhihao Jin, Yingtian Deng and Shiyuan Wu.
