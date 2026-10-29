#!/usr/bin/env python3
"""Post-analysis of saved LFC metrics; no simulations or source-result changes."""
from pathlib import Path
from collections import defaultdict
import csv
import hashlib
import json
import math

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
GRID = ROOT / "results/design_sensitivity"


def read_csv(path):
    rows = []
    with path.open(newline="") as stream:
        for raw in csv.DictReader(stream):
            row = {}
            for key, value in raw.items():
                try:
                    number = float(value)
                    row[key] = number if math.isfinite(number) else None
                except (ValueError, TypeError):
                    row[key] = value
            rows.append(row)
    return rows


def write_csv(name, rows):
    with (OUT / name).open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def passed(row, limits):
    if row["Stable"] != 1:
        return False
    for key, limit in limits.items():
        value = row[key]
        if value is None:
            return False
        if key in ("PM_deg", "GM_dB"):
            if value < limit:
                return False
        elif value > limit:
            return False
    return True


def main():
    original = json.loads((ROOT / "results/research_summary.json").read_text())
    spec = original["specifications"]
    baseline = {
        "Recovery_s": spec["recoveryTime_s"],
        "MaxAbs_Hz": spec["maxAbsoluteDeviationHz"],
        "PM_deg": spec["minPhaseMarginDeg"],
        "GM_dB": spec["minGainMarginDb"],
        "PeakCommand_pu": spec["maxSecondaryCommandPu"],
    }
    register = read_csv(GRID / "design_register.csv")
    rows = read_csv(GRID / "all_plant_cases.csv")
    nominal = {r["DesignID"]: r for r in read_csv(GRID / "nominal_metrics.csv")}
    saved_summary = {r["DesignID"]: r for r in read_csv(GRID / "design_summary.csv")}
    groups = defaultdict(list)
    for row in rows:
        groups[row["DesignID"]].append(row)
    assert len(register) == 36 and len(groups) == 26 and len(rows) == 1118
    for design_id, group in groups.items():
        assert len(group) == 43 and {r["CaseID"] for r in group} == set(range(1, 44))
        assert all(passed(r, baseline) == bool(r["AllTargetsPass"]) for r in group)
        assert sum(passed(r, baseline) for r in group) == saved_summary[design_id]["AllTargetsPassCases"]
    assert max(abs(r["SteadyFrequency_Hz"]) for r in rows) < 1e-10
    # The mathematical zero-frequency equilibrium follows from stable integral
    # control. This residual check is numerical QA, not an extra design threshold.

    # Retrospective policy, explicitly not predeclared before the original study:
    # maximise deterministic all-target case count, then minimise nominal g100.
    # Remaining ties use worst recovery, nominal recovery, and ID, in that order.
    # This introduces priorities but no invented noise acceptance budget.
    def ranking_key(design_id, count):
        recovery = [r["Recovery_s"] for r in groups[design_id]]
        worst = math.inf if None in recovery else max(recovery)
        return (-count, nominal[design_id]["CommandNoiseGain100"], worst,
                nominal[design_id]["Recovery_s"], design_id)

    variations = {
        "Recovery_s": [6, 8, 10],
        "MaxAbs_Hz": [0.020, 0.025, 0.030],
        "PM_deg": [40, 45, 50],
        "GM_dB": [3, 6, 9],
        "PeakCommand_pu": [0.015, 0.020, 0.025],
    }
    scenarios = [("Baseline", "None", None, baseline)]
    for metric, values in variations.items():
        for value in values:
            if value != baseline[metric]:
                scenarios.append((f"{metric}={value:g}", metric, value,
                                  dict(baseline, **{metric: value})))
    scenario_rows, counts_rows, all_scenarios = [], [], []
    baseline_chosen = None
    focus_ids = ["PI_w24_pm65", "PIDF_w24_pm65_N10", "PIDF_w24_pm75_N10"]
    for name, metric, value, limits in scenarios:
        counts = {key: sum(passed(r, limits) for r in group) for key, group in groups.items()}
        ranked = sorted(counts, key=lambda key: ranking_key(key, counts[key]))
        chosen = ranked[0]
        if baseline_chosen is None:
            baseline_chosen = chosen
        maximum = counts[chosen]
        leaders = [key for key in ranked if counts[key] == maximum]
        summary = {
            "Scenario": name, "ChangedMetric": metric, "ChangedValue": value,
            "MaximumPassCases": maximum, "CaseDenominator": 43,
            "PositiveDesignDenominator": 26, "ChosenDesignID": chosen,
            "ChosenNominalNoise100_pu_per_Hz": nominal[chosen]["CommandNoiseGain100"],
            "ChosenNominalRecovery_s": nominal[chosen]["Recovery_s"],
            "ChosenWorstRecovery_s": saved_summary[chosen]["WorstRecovery_s"],
            "SwitchedFromBaseline": chosen != baseline_chosen,
            "MaxCountCandidateIDs": ";".join(leaders),
            "OriginalPI_PassCases": counts[focus_ids[0]],
            "OriginalPIDF_PassCases": counts[focus_ids[1]],
            "BaselineFollowUp_PassCases": counts[focus_ids[2]],
        }
        scenario_rows.append(summary)
        all_scenarios.append(dict(summary, thresholds=limits, counts=counts))
        for reg in register:
            design_id = reg["DesignID"]
            counts_rows.append({
                "Scenario": name, "DesignID": design_id,
                "Status": reg["Status"], "Evaluated": design_id in groups,
                "PassCases": counts.get(design_id),
                "CaseDenominator": 43 if design_id in groups else None,
                "SelectedByRetrospectivePolicy": design_id == chosen,
            })

    leaders = all_scenarios[0]["MaxCountCandidateIDs"].split(";")
    comparisons = []
    failure_details = {}
    for design_id in leaders:
        n, s = nominal[design_id], saved_summary[design_id]
        comparisons.append({
            "DesignID": design_id, "Wc_rad_s": n["Wc_rad_s"],
            "TargetPM_deg": n["TargetPM_deg"], "FilterN": n["FilterN"],
            "PassCases": s["AllTargetsPassCases"], "CaseDenominator": 43,
            "NominalNoise100_pu_per_Hz": n["CommandNoiseGain100"],
            "NominalRecovery_s": n["Recovery_s"],
            "WorstRecovery_s": s["WorstRecovery_s"],
            "WorstPM_deg": s["WorstPM_deg"], "WorstGM_dB": s["WorstGM_dB"],
            "WorstMaxAbs_Hz": s["WorstMaxAbs_Hz"],
            "SelectedByRetrospectivePolicy": design_id == baseline_chosen,
        })
        failure_details[design_id] = [r for r in groups[design_id] if not passed(r, baseline)]

    droop = next(r for r in original["nominal"] if r["name"] == "Droop")
    droop_peak, droop_offset = droop["maxAbsHz"], abs(droop["steadyFrequencyHz"])
    load = original["p"]["loadStep"]
    scales = {
        "status": "Retrospective arithmetic scale checks, not the original derivation, standards or equipment limits.",
        "droopMaximumAbs_Hz": droop_peak, "droopAbsSteadyOffset_Hz": droop_offset,
        "deviationLimit_Hz": baseline["MaxAbs_Hz"],
        "deviationLimitOverDroopPeak": baseline["MaxAbs_Hz"] / droop_peak,
        "requiredReductionFromDroopPeak_percent": 100 * (1 - baseline["MaxAbs_Hz"] / droop_peak),
        "deviationLimitOverDroopOffset": baseline["MaxAbs_Hz"] / droop_offset,
        "recoveryBand_Hz": spec["recoveryBandHz"],
        "bandOverDroopOffset": spec["recoveryBandHz"] / droop_offset,
        "residualOffsetReduction_percent": 100 * (1 - spec["recoveryBandHz"] / droop_offset),
        "commandLimitOverLoadStep": baseline["PeakCommand_pu"] / load,
        "recoveryDeadline_s": baseline["Recovery_s"],
        "recoveryDeadlineInterpretation": "Droop never restores the nominal-frequency band; no finite droop recovery defines the8s deadline.",
        "originalPI_DeadlineSlack_s": baseline["Recovery_s"] - nominal[focus_ids[0]]["Recovery_s"],
        "originalPIDF_DeadlineSlack_s": baseline["Recovery_s"] - nominal[focus_ids[1]]["Recovery_s"],
        "marginLimitsInterpretation": "The45degree/6dB values remain declared comparison criteria, not derived physical limits.",
    }
    source_paths = [GRID / name for name in ("design_register.csv", "all_plant_cases.csv", "nominal_metrics.csv", "design_summary.csv")]
    source_paths += [ROOT / "results/research_summary.json"]
    summary = {
        "scope": "Retrospective saved-metric analysis only; no simulation rerun, gain changes, or original-result edits.",
        "registeredDesigns": 36, "positiveEvaluatedDesigns": 26,
        "unevaluatedInfeasibleDesigns": 10, "casesPerEvaluatedDesign": 43,
        "baselineCriteria": baseline, "recoveryBandUnchanged_Hz": spec["recoveryBandHz"],
        "metricProvenance": "Design-grid margins enumerate polynomial crossings; historical margin approximations are not substituted. Original trajectories and recovery definitions remain unchanged.",
        "baselineMaximumPassCases": max(r["PassCases"] for r in comparisons),
        "baselineFullyQualifyingDesigns": 0, "baselineMaximumCandidates": comparisons,
        "baselineMaximumCandidateFailures": failure_details,
        "policy": "Retrospectively maximise all-target case count, then minimise absolute nominal |CcS(j100)|. Further ties use lower worst recovery, nominal recovery, then DesignID. This is an explicit prioritisation, not an optimisation or noise acceptance budget.",
        "policyChosenDesignID": baseline_chosen,
        "OATvalues": variations, "scenarios": all_scenarios,
        "changedChoices": [s for s in all_scenarios if s["SwitchedFromBaseline"]],
        "unperturbedOtherCriteria": "Each scenario changes exactly one threshold. Stability, fixed +/-0.001Hz recovery band, trajectories, all other thresholds,43-case membership and controller gains remain fixed.",
        "missingRecoveryPolicy": "The one censored original30s row remains a recovery failure for every finite tested deadline; its later diagnostic is not substituted.",
        "thresholdDependenceConclusion": "The baseline candidate persists in8of10 individual perturbations. Tightening maximum deviation to0.020Hz or minimum PM to50degrees selects PIDF_w30_pm75_N10. Loosening PM to40degrees lets four designs pass43/43.",
        "scaleChecks": scales,
        "boundaries": [
            "Selection is among26 fixed positive-gain designs on43 deterministic cases, not a continuous search or probability estimate.",
            "The10 infeasible registered designs remain explicitly unevaluated; they are not assigned fabricated zero-case scores.",
            "The cross-design decision is retrospective and does not replace matched PI/PIDF comparison claims.",
            "Noise100 is one-frequency secondary-measurement-to-command magnitude, not integrated noise variance or a real sensor requirement.",
            "OAT thresholds change classification, not controller trajectories or the recovery band. They do not certify joint threshold variation or physical validity.",
            "Ties and rank switches are disclosed. Candidate preference requires a stated engineering priority and real constraints before deployment.",
        ],
        "validation": {"baselineFlagsAndCountsReproduced": True, "allGroupsHave43UniqueCases": True,
                       "originalDataChanged": False, "simulationsRerun": False},
        "sourcePathBase": "Project root, two levels above results/design_decision_review; all source paths below are relative to that root.",
        "sources": [{"path": path.relative_to(ROOT).as_posix(), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()} for path in source_paths],
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n")
    write_csv("maximum_count_candidates.csv", comparisons)
    write_csv("threshold_scenario_summary.csv", scenario_rows)
    write_csv("all_design_threshold_counts.csv", counts_rows)
    lines = ["RETROSPECTIVE DESIGN DECISION EVIDENCE", "", summary["scope"],
             "36 candidates were registered;26 positive-gain designs were evaluated on43 cases each.",
             "No design passes all 43 cases under inherited criteria. Four tie at 42/43.", "", "Explicit retrospective decision policy:", summary["policy"], "",
             "Maximum-count candidates (all42/43):",
             "Design | g100 (pu/Hz) | nominal recovery(s) | worst recovery(s) | minimum PM(deg)"]
    for c in comparisons:
        lines.append(f"{c['DesignID']} | {c['NominalNoise100_pu_per_Hz']:.6f} | {c['NominalRecovery_s']:.3f} | {c['WorstRecovery_s']:.3f} | {c['WorstPM_deg']:.6f}")
    lines += ["", "The policy selects PIDF(2.4rad/s,75deg,N10), the lowest-noise member of the maximum-count group.",
              "Its sole miss is case31 PM=44.15826768deg,0.84173232deg below45deg. It is not fully qualified.",
              "Faster alternatives exist, so the chosen design is not best in every metric.", "", "One-at-a-time threshold analysis:"]
    for s in scenario_rows:
        lines.append(f"{s['Scenario']}: selected {s['ChosenDesignID']}; maximum {s['MaximumPassCases']}/43; baseline follow-up {s['BaselineFollowUp_PassCases']}/43; original PI/PIDF {s['OriginalPI_PassCases']}/{s['OriginalPIDF_PassCases']}.")
    lines += ["", summary["thresholdDependenceConclusion"],
              "The two switches choose PIDF(3.0rad/s,75deg,N10), whose nominal g100=4.445846pu/Hz, recovery=2.127s and worst recovery=5.372s.",
              "The0.020Hz scenario gives this alternative42/43 versus35/43 for the baseline follow-up. At50deg PM, counts are41/43 versus37/43.",
              "", "Retrospective scale checks; not original derivations or external standards:",
              f"The0.025Hz excursion criterion is {100*scales['deviationLimitOverDroopPeak']:.2f}% of the droop peak0.030697Hz, requiring {scales['requiredReductionFromDroopPeak_percent']:.2f}% reduction.",
              f"The0.001Hz recovery band is {100*scales['bandOverDroopOffset']:.2f}% of the droop offset0.023529Hz, requiring {scales['residualOffsetReduction_percent']:.2f}% offset suppression.",
              "The0.020pu command criterion is twice the0.010pu load. It is not an identified actuator rating.",
              "Droop never returns to the zero-centred band, so it supplies no finite baseline from which8s can be derived.",
              "The original PI/PIDF leave3.345/5.009s of deadline slack; these are observed checks, not proof of the threshold's physical adequacy.",
              "", "Suggested compact assessed-report table: the four maximum-count rows above, with counts stated once in the caption.",
              "Suggested decision: Under the stated retrospective count-first, then lower-noise priority, select2.4/75/N10 for follow-up only. None qualifies against all original cases. The choice changes under stricter excursion or phase-margin criteria.",
              "", "Evidence boundaries:"] + summary["boundaries"]
    print(json.dumps({"baselineMaximum": 42, "selected": baseline_chosen,
                      "switchScenarios": [s["Scenario"] for s in summary["changedChoices"]],
                      "registeredRowsInDerivedCounts": len(counts_rows), "sourceResultsUnmodified": True}, indent=2))


if __name__ == "__main__":
    main()
