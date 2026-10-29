from pathlib import Path
import csv, json, math, cmath, hashlib
from collections import defaultdict
src=Path(__file__).resolve().parents[1]
out=Path(__file__).resolve().parent
files=['design_sensitivity/design_register.csv','design_sensitivity/design_summary.csv','design_sensitivity/nominal_metrics.csv','design_sensitivity/all_plant_cases.csv','design_sensitivity/matched_comparisons.csv','design_decision_review/threshold_scenario_summary.csv']
def read(n):
 rows=[]
 for r in csv.DictReader((src/n).open()):
  row={}
  for k,v in r.items():
   try:
    val=float(v);row[k]=val if math.isfinite(val) else None
   except ValueError: row[k]=v
  rows.append(row)
 return rows
reg={r['DesignID']:r for r in read(files[0])}; saved_summary={r['DesignID']:r for r in read(files[1])}; nom={r['DesignID']:r for r in read(files[2])}; rows=read(files[3]); pairs=read(files[4]); scenarios=read(files[5])
groups=defaultdict(list)
for r in rows:groups[r['DesignID']].append(r)
base={'Recovery_s':8,'MaxAbs_Hz':.025,'PeakCommand_pu':.020,'PM_deg':45,'GM_dB':6}
def passed(r,limits):
 return r['Stable']==1 and all(r[k] is not None and (r[k]>=v if k in ['PM_deg','GM_dB'] else r[k]<=v) for k,v in limits.items())
assert len(rows)==1118 and len(groups)==26 and len(reg)==36
for k,g in groups.items():
 assert len(g)==43 and {r['CaseID'] for r in g}==set(range(1,44))
 assert all(passed(r,base)==bool(r['AllTargetsPass']) for r in g)
 assert sum(passed(r,base) for r in g)==saved_summary[k]['AllTargetsPassCases']
# Check stored pair metrics by joining their two nominal records.
eligible=[]
for p in pairs:
 if p['Eligible']!=1:continue
 pi,pid=nom[p['PI_DesignID']],nom[p['PIDF_DesignID']]
 assert reg[p['PI_DesignID']]['MatchedEligible']==1 and reg[p['PIDF_DesignID']]['MatchedEligible']==1
 for pk,nk in [('NominalRecovery_s','Recovery_s'),('NominalIAE_Hz_s','IAE_Hz_s'),('Noise100','CommandNoiseGain100')]:
  assert math.isclose(p['PI_'+pk],pi[nk],abs_tol=1e-12)
  assert math.isclose(p['PIDF_'+pk],pid[nk],abs_tol=1e-12)
 eligible.append({'PI_DesignID':p['PI_DesignID'],'PIDF_DesignID':p['PIDF_DesignID'],'wc_rad_s':p['Wc_rad_s'],'PM_deg':p['TargetPM_deg'],'N':p['FilterN'],'recovery_advantage_s':pi['Recovery_s']-pid['Recovery_s'],'IAE_reduction_percent':100*(1-pid['IAE_Hz_s']/pi['IAE_Hz_s']),'noise_ratio_PIDF_over_PI':pid['CommandNoiseGain100']/pi['CommandNoiseGain100']})
assert len(eligible)==15
matched_summary={}
for key in ['recovery_advantage_s','IAE_reduction_percent','noise_ratio_PIDF_over_PI']:
 lo=min(eligible,key=lambda r:r[key]);hi=max(eligible,key=lambda r:r[key])
 matched_summary[key]={'minimum':lo[key],'maximum':hi[key],'minimum_design':lo['PIDF_DesignID'],'maximum_design':hi['PIDF_DesignID']}
grouppairs=defaultdict(list)
for p in eligible:grouppairs[(p['wc_rad_s'],p['PM_deg'])].append(p)
matched_compact=[]
for (wc,pm),g in grouppairs.items():
 assert len(g)==3
 matched_compact.append({'wc_rad_s':wc,'PM_deg':pm,'N_values':[r['N'] for r in g],**{key:[min(r[key] for r in g),max(r[key] for r in g)] for key in ['recovery_advantage_s','IAE_reduction_percent','noise_ratio_PIDF_over_PI']}})
# Reclassify all 11 threshold scenarios from original saved per-case metrics.
thresholds=[]
for sc in scenarios:
 lim=dict(base)
 if sc['ChangedMetric']!='None':lim[sc['ChangedMetric']]=sc['ChangedValue']
 counts={k:sum(passed(r,lim) for r in g) for k,g in groups.items()}
 def rank(k):
  vals=[r['Recovery_s'] for r in groups[k]]
  return (-counts[k],nom[k]['CommandNoiseGain100'],math.inf if None in vals else max(vals),nom[k]['Recovery_s'],k)
 ranking=sorted(groups,key=rank);chosen=ranking[0]
 assert chosen==sc['ChosenDesignID']
 assert counts[chosen]==sc['MaximumPassCases']
 assert counts['PIDF_w24_pm75_N10']==sc['BaselineFollowUp_PassCases']
 assert counts['PI_w24_pm65']==sc['OriginalPI_PassCases']
 assert counts['PIDF_w24_pm65_N10']==sc['OriginalPIDF_PassCases']
 thresholds.append({'scenario':sc['Scenario'],'thresholds':lim,'selected':chosen,'selected_pass_count':counts[chosen],'candidate_A_pass_count':counts['PIDF_w24_pm75_N10'],'original_PI_pass_count':counts['PI_w24_pm65'],'original_PIDF_pass_count':counts['PIDF_w24_pm65_N10'],'maximum_count_ties':[k for k in ranking if counts[k]==counts[chosen]],'denominator':43})
ids=['PI_w24_pm65','PIDF_w24_pm65_N10','PIDF_w24_pm75_N10']
comparison=[]
for k in ids:
 g=groups[k]
 worst={}
 for metric in ['Recovery_s','MaxAbs_Hz','PeakCommand_pu','PM_deg','GM_dB','Ms','CommandNoiseGain100','IAE_Hz_s']:
  row=(min if metric in ['PM_deg','GM_dB'] else max)(g,key=lambda r:r[metric])
  worst[metric]={'value':row[metric],'case':int(row['CaseID'])}
 comparison.append({'DesignID':k,'gains':{q:reg[k][q] for q in ['Kp','Ki','Kd','Tf']},'nominal':{q:nom[k][q] for q in ['Recovery_s','MaxAbs_Hz','IAE_Hz_s','PeakCommand_pu','PM_deg','GM_dB','Ms','CommandNoiseGain100']},'worst_43':worst,'pass_count':sum(passed(r,base) for r in g),'stable_count':sum(r['Stable']==1 for r in g),'failures':[{q:r[q] for q in ['CaseID','Tg_factor','Tt_factor','Tp_factor','Ksys_factor','R_factor','Recovery_s','MaxAbs_Hz','PeakCommand_pu','PM_deg','GM_dB']} for r in g if not passed(r,base)]})
a=reg[ids[-1]];wc=a['Wc_rad_s'];pm=a['TargetPM_deg'];N=a['FilterN'];s=1j*wc
P=120/(.48*s**3+7.624*s**2+20.38*s+51)
Cstar=cmath.exp(1j*math.radians(pm-180))/P
Tf=1/(N*wc);ad=wc*wc*Tf/(1+(wc*Tf)**2);bd=wc/(1+(wc*Tf)**2)
Kd=(Cstar.real+Cstar.imag)/(ad+bd);Kp=Cstar.real-ad*Kd;Ki=wc*Kp
independent={'Kp':Kp,'Ki':Ki,'Kd':Kd,'Tf':Tf}
assert max(abs(a[k]-v) for k,v in independent.items())<1e-13
C=a['Kp']+a['Ki']/s+a['Kd']*s/(1+a['Tf']*s)
units={'Kp':'pu/Hz','Ki':'pu/(Hz s)','Kd':'pu s/Hz','Tf':'s','g100':'pu/Hz'}
output={'scope':'Read-only saved CSV aggregation and algebra check; no simulations or source changes. All precision reflects saved numerical output, not measurement accuracy.','inputs':[{'path':'results/'+f,'sha256':hashlib.sha256((src/f).read_bytes()).hexdigest()} for f in files],'validation':{'rows':len(rows),'designs_registered':len(reg),'designs_evaluated':len(groups),'cases_per_design':43,'baseline_flags_and_counts_match':True,'threshold_summary_independently_matches':True,'eligible_pair_nominal_joins_match':True,'final_gains_independent_algebra':independent,'final_point_match_residual':abs(C*P-cmath.exp(1j*math.radians(pm-180)))},'gain_units':units,'comparison':comparison,'matched_pair_summary':matched_summary,'matched_compact_5rows':matched_compact,'all15_pairs':eligible,'threshold_scenarios':thresholds,'same_A_choice_in_10_changes':sum(r['selected']==ids[-1] for r in thresholds[1:]),'candidate_noise_increase_percent_vs_original_PIDF':100*(nom[ids[-1]]['CommandNoiseGain100']/nom[ids[1]]['CommandNoiseGain100']-1),'candidate_case31_PM_shortfall_deg':45-comparison[-1]['worst_43']['PM_deg']['value']}
(out/'evidence.json').write_text(json.dumps(output,indent=2,allow_nan=False)+'\n')
lines=['V4 COMPACT EVIDENCE — independent saved-metric audit','',output['scope'],'','All source paths and SHA256 hashes are in evidence.json.','Final development candidate: PIDF, wc=2.4 rad/s, target PM=75 deg, N=10.','Gains (rounded for report): Kp=0.30894426 pu/Hz; Ki=0.74146622 pu/(Hz s); Kd=0.06763856 pu s/Hz; Tf=0.04166667 s.','43/43 stable; 42/43 pass all declared thresholds. Case31 fails only PM: 44.15826768 deg, shortfall0.84173232 deg.','', 'RECOMMENDED TABLE 1: original PI65 / original PIDF65 / development PIDF75', 'Metric | PI65 | PIDF65 | PIDF75']
for key,fmt in [('Recovery_s','.3f'),('IAE_Hz_s','.6f'),('MaxAbs_Hz','.6f'),('PeakCommand_pu','.6f'),('PM_deg','.3f'),('GM_dB','.3f'),('Ms','.3f'),('CommandNoiseGain100','.3f')]:
 lines.append('Nominal '+key+' | '+' | '.join(format(r['nominal'][key],fmt) for r in comparison))
for key,fmt in [('Recovery_s','.3f'),('MaxAbs_Hz','.6f'),('PeakCommand_pu','.6f'),('PM_deg','.3f'),('GM_dB','.3f')]:
 lines.append('Worst43 '+key+' | '+' | '.join(format(r['worst_43'][key]['value'],fmt) for r in comparison))
lines+=['All-target passes | '+' | '.join(str(r['pass_count'])+'/43' for r in comparison),'','Suggested compact version: gains in candidate caption/box; table rows nominal/worst recovery, nominalIAE, worstexcursion, minimumPM, nominalg100, pass count. State per-row extremes need not occur in same case.','', 'RECOMMENDED TABLE 2: 15 nominal matched pairs grouped by shared(wc,PM); N=5,10,20 in each row.', 'wc | PM | recovery advantage range(s) | IAE reduction range(%) | g100 ratio range(PIDF/PI)']
for r in matched_compact:
 vals=['–'.join(f'{v:.3f}' for v in r[k]) for k in ['recovery_advantage_s','IAE_reduction_percent','noise_ratio_PIDF_over_PI']]
 lines.append(f"{r['wc_rad_s']:g} | {r['PM_deg']:g} | "+' | '.join(vals))
lines+=['Across all15: nominal recovery improves0.007–16.616s; nominalIAE reduction0.884–85.358%; noise ratio1.069–15.795. These are observed bounds, not confidence intervals or independent statistical trials.','The largest benefit includes PI(3.0,65) with18.591s nominal recovery, which fails the8s deadline. Do not present that extreme as a typical gain.','', 'RECOMMENDED TABLE 3: threshold sensitivity, one threshold changed at a time.', 'Scenario | selected A/B | selected pass count | A pass count | original PI/PIDF counts']
for r in thresholds:
 short='A' if r['selected']==ids[-1] else 'B'
 lines.append(f"{r['scenario']} | {short} | {r['selected_pass_count']}/43 | {r['candidate_A_pass_count']}/43 | {r['original_PI_pass_count']}/{r['original_PIDF_pass_count']}")
lines+=['A=PIDF(2.4,75,N10); B=PIDF(3.0,75,N10). Other thresholds, recovery band and trajectories remain unchanged.','For a5-row table, combine paired threshold values, and show corresponding selected A/B and counts in same order.','A remains selected in8of10 changes. It is not joint-threshold, continuous-parameter or physical robustness.','Original65-degree anti-windup/saturation results do not validate A.','', 'CANDIDATE EXACT FAILURE ROW',json.dumps(comparison[-1]['failures'][0],indent=2),'','Checks: all1118 baseline flags/counts, all11 threshold scenarios, all15 pair joins and selected gain algebra verified. No input files modified.']
(out/'evidence.txt').write_text('\n'.join(lines)+'\n')
print(json.dumps({'output':str(out/'evidence.json'),'candidate':comparison[-1],'summary':matched_summary,'grouped':matched_compact,'samechoice':output['same_A_choice_in_10_changes'],'algebra':output['validation']['final_point_match_residual']},indent=2))
