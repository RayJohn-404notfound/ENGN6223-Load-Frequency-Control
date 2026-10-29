# Results and Interpretation

These results describe a small-deviation, single-area, non-reheat benchmark with primary droop retained inside the secondary-loop plant. They are simulation and saved-data results, not measurements from a physical generator.

## Model and tests

The three states are governor output, mechanical-power deviation and frequency deviation. The parameters are adopted from Anwar and Pan's Example 1:

| Parameter | Value |
|---|---:|
| Governor time constant | 0.08 s |
| Turbine time constant | 0.30 s |
| System time constant | 20 s |
| System gain | 120 Hz/pu |
| Droop | 2.4 Hz/pu |

The nominal test applies a 0.01 pu load increase at 1 s, starting from zero deviations. Recovery is measured after the disturbance: frequency must enter the ±0.001 Hz band and remain there through the 30 s simulation endpoint.

The declared comparison criteria are recovery ≤8 s, maximum frequency deviation ≤0.025 Hz, peak secondary command ≤0.020 pu, phase margin ≥45° and gain margin ≥6 dB, with stable dynamics and zero steady frequency error. These are study choices, not grid-code requirements or identified equipment limits.

The 43 deterministic plant cases comprise the nominal model, 32 parameter corners and 10 individual endpoints. Parameters vary by ±20%, with each controller's gains held fixed. Case counts are not probabilities or a proof over a continuous uncertainty region.

## Original comparison and development candidate

Both original controllers target 2.4 rad/s crossover and 65° phase margin. PIDF uses a derivative filter with its pole at ten times crossover. The later candidate targets 75°, so it is a separate development choice.

| Metric | Original PI | Original PIDF | Selected PIDF candidate |
|---|---:|---:|---:|
| Nominal phase target | 65° | 65° | 75° |
| Nominal recovery (s) | 4.655 | 2.991 | 2.788 |
| Worst recovery across 43 cases (s) | 13.788 | 6.706 | 4.376 |
| Nominal command-noise gain at 100 rad/s (pu/Hz) | 0.294 | 1.100 | 1.879 |
| Minimum phase margin across 43 cases | 20.299° | 32.181° | 44.158° |
| Cases passing every criterion | 34/43 | 36/43 | 42/43 |

The selected PIDF has `Kp = 0.30894426`, `Ki = 0.74146622`, `Kd = 0.06763856` and `Tf = 0.04166667 s`, in parallel form:

```text
C(s) = Kp + Ki/s + Kd*s/(1 + Tf*s)
```

With frequency error in Hz and command in pu, the gain units are pu/Hz, pu/(Hz·s) and pu·s/Hz, respectively.

The grid registers 36 designs; 26 have positive required gains and are evaluated. Ten exclusions remain recorded rather than receiving invented scores. Four designs tie at 42/43. A retrospective rule first maximises the passing-case count, then minimises nominal noise gain at 100 rad/s among tied designs. It selects the candidate above, which still fails case 31 phase margin by 0.842°.

The choice persists under eight of ten individual threshold changes. Stricter frequency-excursion or phase-margin limits select another candidate. The rule is a stated priority, not global optimisation or a physical noise budget.

## Why additional analyses matter

- **Recovery exceptions:** Original PIDF recovers later in cases 15, 30 and 31. In case 30, weaker cancellation between response modes produces an extra recovery-band exit despite faster dominant decay. Increasing the band from 1.0 to 1.1 mHz reverses that ranking.
- **Command saturation:** For the original 65° PIDF, anti-windup reduces recovery after load removal from 17.485 to 4.755 s. The command bound stays at 0.015 pu. The tracking time is a fixed 0.5 s; no tuning sweep establishes its optimality.
- **Measurement noise:** At 100 rad/s, the original PIDF has about 3.74 times PI's secondary-noise-to-command gain. This is a transfer-function calculation at one frequency, not measured sensor noise or total noise variance.

## Evidence map

| Question | Main saved evidence |
|---|---|
| What happens in the original comparison? | `results/controller_gains.csv`, `nominal_metrics.csv`, `parameter_sensitivity.csv` |
| Why do some recovery rankings reverse? | `results/deepening/` |
| How are margin calculations checked? | `results/review_revision/` |
| What happens with command limits? | `results/saturation_metrics.csv` |
| Which designs were evaluated or excluded? | `results/design_sensitivity/design_register.csv`, `all_plant_cases.csv`, `design_summary.csv` |
| Why was the candidate selected? | `results/design_decision_review/maximum_count_candidates.csv`, `threshold_scenario_summary.csv` |
| How are compact summaries and figures regenerated? | Python scripts in `results/report_evidence_summary/` and `results/report_compact/` |

## Limits and next steps

Numerical agreement between formulations checks implementation consistency. It does not validate the benchmark against an actual plant. The model omits governor deadband, generation-rate limits, valve travel, delays and network dynamics.

Noise enters only the secondary measurement; primary droop measurement is ideal. Noise and saturation are considered separately. The original anti-windup results do not validate the selected 75° candidate. Its next tests should combine noise and saturation, vary anti-windup tuning, and use identified equipment constraints and parameter correlations.
