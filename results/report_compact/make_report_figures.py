"""Regenerate report figures from existing result tables; no simulations or retuning.
Run with Python, NumPy and Matplotlib. All inputs are relative to this script.
The vector PDFs are designed for a 7.05-inch two-column report width.
"""
from pathlib import Path
import hashlib
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT = Path(__file__).resolve().parent
RESULTS = OUT.parent
PI, PIDF, BASE = '#1764A0', '#C14C18', '#66717A'
plt.rcParams.update({'font.family':'DejaVu Sans', 'font.size':8.5,
    'axes.labelsize':8.5, 'axes.titlesize':9, 'legend.fontsize':8,
    'xtick.labelsize':8, 'ytick.labelsize':8, 'lines.linewidth':1.15,
    'axes.linewidth':0.6, 'grid.linewidth':0.45, 'grid.alpha':0.35,
    'pdf.fonttype':42, 'ps.fonttype':42, 'savefig.facecolor':'white'})
inputs = {}
def read(name):
    p = RESULTS/name
    inputs[name] = hashlib.sha256(p.read_bytes()).hexdigest()
    return np.genfromtxt(p, delimiter=',', names=True, encoding='utf-8')
def style(ax):
    ax.grid(True)
    ax.spines[['top','right']].set_visible(False)
def output(fig, name):
    fig.savefig(OUT/(name+'.pdf'))
    fig.savefig(OUT/(name+'.png'), dpi=200)
    plt.close(fig)
def pair():
    fig, ax=plt.subplots(1,2,figsize=(7.05,2.35),layout='constrained')
    for a in ax: style(a)
    return fig,ax

x=read('nominal_response.csv')
fig,ax=pair()
for name,col,colour,ls in [('Droop','Droop_df_Hz',BASE,'-'),('PI','PI_df_Hz',PI,'--'),('PIDF','PIDF_df_Hz',PIDF,'-')]:
    ax[0].plot(x['Time_s'],x[col],color=colour,ls=ls,label=name)
for y in [-.001,.001]: ax[0].axhline(y,color='#999999',lw=.65,ls=':')
ax[0].set(xlim=(0,15),ylim=(-.035,.005),xlabel='Time (s)',ylabel='Frequency deviation (Hz)',title='(a) Frequency recovery')
ax[0].legend(loc='lower right',frameon=False,ncol=3,handlelength=1.6,columnspacing=.8)
for name,colour,ls in [('PI',PI,'--'),('PIDF',PIDF,'-')]:
    ax[1].plot(x['Time_s'],x[name+'_u_pu'],color=colour,ls=ls,label=name)
ax[1].axhline(.01,color='#777777',ls=':',lw=.7)
ax[1].set(xlim=(0,15),ylim=(0,.012),xlabel='Time (s)',ylabel='Secondary command (pu)',title='(b) Command demand')
ax[1].legend(loc='lower right',frameon=False)
for a in ax: a.axvline(1,color='#999999',ls=':',lw=.7); a.set_xticks([0,5,10,15])
output(fig,'nominal')

b=read('bode_data.csv')
fig,ax=pair()
for name,colour,ls in [('PI',PI,'--'),('PIDF',PIDF,'-')]:
    ax[0].semilogx(b['omega_rad_s'],b[name+'_magnitude_dB'],color=colour,ls=ls,label=name)
    ax[1].semilogx(b['omega_rad_s'],b[name+'_phase_deg'],color=colour,ls=ls,label=name)
for a in ax:
    a.set_xlim(.01,100); a.set_xlabel('Angular frequency (rad/s)')
    a.axvline(2.4,color='#777777',ls=':',lw=.8);a.legend(loc='lower left',frameon=False)
ax[0].axhline(0,color='#888888',ls=':',lw=.7)
ax[0].set(ylim=(-100,50),ylabel='Loop magnitude (dB)',title='(a) Magnitude')
ax[1].axhline(-180,color='#888888',ls=':',lw=.7)
ax[1].set(ylim=(-300,-50),ylabel='Loop phase (degrees)',title='(b) Phase')
output(fig,'bode')

s=read('deepening/sensitivity_peak_data.csv')
p=RESULTS/'deepening/phase_margin_summary.json'
inputs[str(p.relative_to(RESULTS))]=hashlib.sha256(p.read_bytes()).hexdigest()
meta=json.loads(p.read_text())['nominal']
fig,ax=pair()
for i,(name,colour,ls) in enumerate([('PI',PI,'--'),('PIDF',PIDF,'-')]):
    ax[0].semilogx(s['omega_rad_s'],s[name+'_S_magnitude'],color=colour,ls=ls,label=name)
    m=meta[i]
    ax[0].plot(m['Ms_peak_rad_s'],m['Ms_peak'],'o',color=colour,ms=4)
    L=10**(b[name+'_magnitude_dB']/20)*np.exp(1j*np.deg2rad(b[name+'_phase_deg']))
    ax[1].plot(L.real,L.imag,color=colour,ls=ls,label=name)
    ax[1].plot(m['closestLoopReal'],m['closestLoopImag'],'o',color=colour,ms=4)
    ax[1].plot([-1,m['closestLoopReal']],[0,m['closestLoopImag']],color=colour,ls=':',lw=.75)
ax[0].plot(2.4,meta[0]['sensitivityAtCrossover'],'s',color='black',ms=3.5,label='Shared crossover')
ax[0].axvline(2.4,color='#777777',ls=':',lw=.7)
ax[0].set(xlim=(.3,30),ylim=(0,2.1),xlabel='Angular frequency (rad/s)',ylabel=r'$|S(j\omega)|$',title='(a) Sensitivity peak differs')
ax[0].legend(loc='lower right',frameon=False)
ax[1].plot(-1,0,'x',color='#A80000',ms=6,mew=1.5,label='Critical point')
ax[1].plot(meta[0]['commonLoopReal'],meta[0]['commonLoopImag'],'s',color='black',ms=3.5,label='Shared crossover')
ax[1].set(xlim=(-1.2,.2),ylim=(-1.2,.3),xlabel=r'Real $L(j\omega)$',ylabel=r'Imaginary $L(j\omega)$',title='(b) Positive-frequency detail')
ax[1].set_aspect('equal',adjustable='box')
ax[1].text(-1.16,.17,'Critical point -1',fontsize=8,color='#A80000')
output(fig,'sensitivity_nyquist')

fig,ax=pair()
for name,colour in [('PI',PI),('PIDF',PIDF)]:
    for aw in [0,1]:
        d=read(f'saturation_{name}_AW{aw}.csv')
        label=f'{name}, AW '+('on' if aw else 'off')
        for a,key in zip(ax,['FrequencyDeviation_Hz','Command_pu']):
            a.plot(d['Time_s'],d[key],color=colour,ls='-' if aw else '--',label=label)
for a in ax:
    a.set_xlim(0,50);a.set_xlabel('Time (s)')
    for t in [5,20]:a.axvline(t,color='#888888',ls=':',lw=.7)
ax[0].set(ylim=(-.08,.06),ylabel='Frequency deviation (Hz)',title='(a) Recovery after load removal')
ax[1].set(ylim=(-.011,.016),ylabel='Applied command (pu)',title='(b) Limited secondary command')
ax[0].legend(loc='lower right',frameon=False,ncol=2,columnspacing=.7,handlelength=1.5)
ax[1].text(21.5,.0092,'Load off',fontsize=8,color='#555555')
output(fig,'antiwindup')

(OUT/'figure_provenance.json').write_text(json.dumps({'purpose':'Report figures regenerated from saved result tables; no rerun or retuning', 'input_sha256':inputs,'figure_width_inches':7.05,'body_label_font_pt':8.5,'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},indent=2)+'\n')
print('Created four vector PDF figures and PNG previews from unchanged source tables.')
