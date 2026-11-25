#!/usr/bin/env python3
"""Regenerate manuscript tables and figures from the verified full-results archive.
Usage: python make_paper_assets.py results/full audit/derived manuscript
"""
from pathlib import Path
import sys
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
R=Path(sys.argv[1]); A=Path(sys.argv[2]); M=Path(sys.argv[3])
T=M/'tables'; F=M/'figures'; T.mkdir(parents=True,exist_ok=True); F.mkdir(parents=True,exist_ok=True)
def rd(root,name):return pd.read_csv(root/'tables'/name)
U={6:R/'main_no_vix',7:R/'extended_with_vix'}
S=pd.read_csv(A/'score_summary.csv')
P=pd.read_csv(A/'paired_block20_summary.csv')
LAB=dict(gaussian='Gaussian',t_copula='Student-$t$',static_vine='Static vine',regime_vine='Regime vine',dcc_gaussian='Gaussian DCC',rolling_vine='Rolling (tuned)',rolling_matched='Rolling (matched)',WRR_direct='WRR (direct)',WRR_raw='Raw PITs',WRR_no_tail='No tail bonus',WRR_no_recency='No recency',WRR_no_regime='No state bonus',WRR_no_vix_distance='No VIX proximity',WRR_other_family='Other family pool',WRR_legacy_resampled_B1='Raw PITs, $B=1$',WRR_resampled_B1='Rank-adjusted $B=1$',WRR_resampled_B5='Rank-adjusted $B=5$',WRR_resampled_B20='Rank-adjusted $B=20$',WRR_threshold_55_80='Thresholds 55/80',WRR_threshold_65_90='Thresholds 65/90')
base=['WRR_direct','rolling_matched','rolling_vine','dcc_gaussian','t_copula','static_vine','regime_vine','gaussian']
ablation=['WRR_direct','WRR_raw','WRR_no_tail','WRR_no_recency','WRR_no_regime','WRR_no_vix_distance','WRR_other_family','WRR_legacy_resampled_B1','WRR_resampled_B1','WRR_resampled_B5','WRR_resampled_B20','WRR_threshold_55_80','WRR_threshold_65_90']
def fx(x,n=4):
    return f'{x:.{n}f}'
def pfx(x):
    if pd.isna(x): return '--'
    if x<.0001:return '$<0.0001$'
    return f'{x:.4f}'
def tab(name,align,headers,rows):
    s=['\\begin{tabular}{'+align+'}','\\toprule']+headers+['\\midrule']+[' & '.join(map(str,row))+r' \\' for row in rows]+['\\bottomrule','\\end{tabular}']
    (T/(name+'.tex')).write_text('\n'.join(s)+'\n',encoding='utf-8')
def modelval(spec,d,model,metric):return S.loc[(S.spec==spec)&(S.dimension==d)&(S.model==model),metric].iloc[0]
desc=rd(R,'descriptive_returns.csv')
assets=dict(XLE='Conventional energy',ICLN='Clean energy',TAN='Solar energy',XLF='Financial sector',GLD='Gold',UUP='US dollar',VIX='Volatility index change')
tab('descriptive','llrrrrr',[r'Asset & Exposure & Mean & SD & 5th pct. & Median & 95th pct. \\'],[[r.asset,assets[r.asset],*[fx(100*v,3) for v in [r.mean,r.sd,r.q05,r.median,r.q95]]] for r in desc.itertuples()])
cal=pd.read_csv(A/'pit_calibration_checked.csv').query("spec=='primary' and split=='test' and q==0.05")
tab('pit','lrrrr',[r'Asset & $n$ & Hit rate (\%) & 95\% interval (\%) & Holm $p$ \\'],[[r.asset,r.n,fx(100*r.hit_rate,2),f'[{100*r.hit_ci_low:.2f}, {100*r.hit_ci_high:.2f}]',pfx(r.holm_p_documented)] for r in cal.itertuples()])
reg=rd(U[6],'regimes.csv')
rg=pd.crosstab(reg.split,reg.regime)
tab('states','lrrrr',[r'Period & Calm & Medium & Stress & Total \\'],[[s.title(),*rg.loc[s,['calm','medium','stress']].tolist(),rg.loc[s].sum()] for s in ['train','validation','test']])
def compare(name,models):
    heads=[r' & \multicolumn{2}{c}{Log score $\uparrow$} & \multicolumn{2}{c}{Energy $\downarrow$} & \multicolumn{2}{c}{Variogram $\downarrow$} & \multicolumn{2}{c}{Brier $\downarrow$} \\',r'\cmidrule(lr){2-3}\cmidrule(lr){4-5}\cmidrule(lr){6-7}\cmidrule(lr){8-9}',r'Method & 6D & 7D & 6D & 7D & 6D & 7D & 6D & 7D \\']
    rows=[]
    for model in models:
        rr=[LAB[model]]
        for metric,dp,maxim in [('mean_log_score',4,True),('energy_common',4,False),('variogram_common',4,False),('brier',6,False)]:
            for d in [6,7]:
                val=modelval('primary',d,model,metric)
                vals=[modelval('primary',d,m,metric) for m in models]
                s=fx(val,dp)
                if val==(max(vals) if maxim else min(vals)):s=r'\textbf{'+s+'}'
                rr.append(s)
        rows.append(rr)
    tab(name,'lrrrrrrrr',heads,rows)
compare('benchmarks',base);compare('ablations',ablation)
# Paired comparisons: two dimensions, positive direct-minus-comparator loss.
rows=[]
for model,metric in [(m,'nll') for m in ['rolling_matched','rolling_vine','dcc_gaussian','t_copula','WRR_raw','WRR_no_tail','WRR_resampled_B20']]+[(m,'variogram_common') for m in ['WRR_raw','WRR_no_tail']]:
    rr=[LAB[model],'NLL' if metric=='nll' else 'Variogram']
    for d in [6,7]:
        r=P.loc[(P.spec=='primary')&(P.dimension==d)&(P.metric==metric)&(P.benchmark==model)].iloc[0]
        rr+=[fx(r.difference,4),f'[{r.ci_low:.4f}, {r.ci_high:.4f}]',pfx(r.p_holm)]
    rows.append(rr)
tab('inference','llrcrrcr',[r' & & \multicolumn{3}{c}{6D} & \multicolumn{3}{c}{7D} \\',r'\cmidrule(lr){3-5}\cmidrule(lr){6-8}',r'Comparator & Loss & Diff. & 95\% interval & $p_H$ & Diff. & 95\% interval & $p_H$ \\'],rows)
# Same-target extended-minus-main comparisons.
cross=rd(R,'common_six_cross_universe.csv')
rn={'energy_common':'Energy','variogram_common':'Variogram','brier':'Brier','quantile_loss':'Quantile loss'}
tab('cross','lrrr',[r'Loss & Difference & 95\% interval & Holm $p$ \\'],[[rn[r.metric],fx(r.extended_minus_main,6),f'[{r.ci_low:.6f}, {r.ci_high:.6f}]',pfx(r.p_holm)] for r in cross.query("model=='WRR_direct'").itertuples()])
# Regime heterogeneity (same-day log score differences).
rg=pd.read_csv(A/'regime_logscore_contrasts.csv'); rsp={'primary':'GARCH, untrimmed','empirical_margin':'Empirical margins','train_winsor':'Training bounds'}
rows=[]
for spec in rsp:
    rr=[rsp[spec]]
    for d in [6,7]:
        for s in ['calm','medium','stress']:
            v=rg.query('spec==@spec and dimension==@d and regime==@s').WRR_minus_matched_log.iloc[0]
            rr.append(f'{v:+.4f}')
    rows.append(rr)
tab('regime_contrasts','lrrrrrr',[r'& \multicolumn{3}{c}{6D} & \multicolumn{3}{c}{7D} \\',r'\cmidrule(lr){2-4}\cmidrule(lr){5-7}',r'Margins / preprocessing & Calm & Medium & Stress & Calm & Medium & Stress \\'],rows)
# Rank diagnostics and numerical uncertainty.
unc=pd.read_csv(A/'uncertainty_summary.csv')
rows=[]
for d in [6,7]:
    diag=rd(U[d],'refit_weight_diagnostics.csv'); dd=diag.query("model=='WRR_direct'"); no=diag.query("model=='WRR_no_tail'"); ens=diag.query("model=='WRR_resampled_B20'"); un=unc.query('dimension==@d')
    cs=rd(U[d],'covar_mc_summary.csv')
    rows.append([d,int(dd.n_window.iloc[0]),fx(dd.ess.mean(),1),fx((dd.ess/dd.n_window).min()*100,2),fx(dd.weighted_ks_before.mean()),fx(dd.weighted_ks_after.mean()),fx(no.ess.mean(),1),fx(ens.unique_mean.mean(),1)])
tab('weights','rrrrrrrr',[r'$d$ & $W$ & Mean ESS & Min. ESS (\%) & Before KS & After KS & No-tail ESS & $B=20$ unique \\'],rows)
rows=[]
for d in [6,7]:
    a=unc.query('dimension==@d');cs=rd(U[d],'covar_mc_summary.csv')
    rows.append([d,1791,fx(100*(a.boot_high-a.boot_low).mean(),3),fx(100*a.boot_sd.mean(),3),*[fx(100*a['mc_sd_'+str(n)].mean(),3) for n in [2000,5000,20000]],int(cs.min_joint_tail.min())])
tab('uncertainty','rrrrrrrr',[r' & & \multicolumn{2}{c}{Local refits} & \multicolumn{3}{c}{MC SD (probability, pp)} & \\',r'\cmidrule(lr){3-4}\cmidrule(lr){5-7}',r'$d$ & Successes & Width (pp) & SD (pp) & 2,000 & 5,000 & 20,000 & Min. tail $n$ \\'],rows)
# Robustness and selection.
rows=[]
for spec in rsp:
    for d in [6,7]:
        rr=[rsp[spec],d]+[fx(modelval(spec,d,m,'mean_log_score'),4) for m in ['WRR_direct','rolling_matched','rolling_vine','dcc_gaussian']]
        rr.append(fx(modelval(spec,d,'WRR_direct','mean_log_score')-modelval(spec,d,'rolling_matched','mean_log_score'),4))
        rr.append(fx(modelval(spec,d,'WRR_direct','probability')*100,2))
        rr.append(fx(modelval(spec,d,'WRR_direct','observed')*100,2))
        rows.append(rr)
tab('robustness','lrrrrrrrr',[r'Specification & $d$ & WRR & Matched & Tuned & DCC & WRR$-$matched & $\bar p$ (\%) & Obs. (\%) \\'],rows)
se=pd.read_csv(A/'validation_selection.csv')
tab('selection','lrlrrlr',[r'Specification & $d$ & Model & $W$ & Refit & Pool & Val. log score \\'],[[rsp[r.spec],r.dimension,LAB[r.model],r.window,r.refit,r.family,fx(r.validation_log_score)] for r in se.itertuples()])
# Additional tail event thresholds.
gr=pd.read_csv(A/'tail_grid_summary.csv').query("spec=='primary'")
rows=[]
for (q,k),g in gr.groupby(['q','k']):
    row=[fx(q,2),k,int(round(502*g.observed.iloc[0]))]
    for d in [6,7]:
        for m in ['WRR_direct','rolling_matched']:
            row.append(fx(g.loc[(g.dimension==d)&(g.model==m),'brier'].iloc[0],6))
    rows.append(row)
tab('tailgrid','rrrcccc',[r'& & & \multicolumn{2}{c}{6D Brier} & \multicolumn{2}{c}{7D Brier} \\',r'\cmidrule(lr){4-5}\cmidrule(lr){6-7}',r'$q$ & $k$ & Events & WRR & Matched & WRR & Matched \\'],rows)
# Risk summaries: same return target throughout.
rows=[]
for d in [6,7]:
    for m in ['WRR_direct','WRR_no_tail','rolling_matched','rolling_vine','dcc_gaussian','t_copula']:
        vals=[modelval('primary',d,m,k) for k in ['violation','quantile_loss','es_identification']]
        rows.append([d,LAB[m],int(round(502*vals[0])),fx(vals[0]*100,2),fx(vals[1]*10000,3),fx(vals[2]*10000,3)])
tab('portfolio','rlrrrr',[r'$d$ & Model & Violations & Rate (\%) & Quantile loss (bp) & ES moment (bp) \\'],rows)
# Observed finite-threshold network, selected pairs.
net=rd(U[6],'finite_threshold_network.csv')
rows=[]
for i,j in [('XLE','XLF'),('ICLN','TAN'),('ICLN','XLF'),('TAN','XLF'),('XLE','GLD')]:
    rr=[i+' $\\to$ '+j]
    for s in ['calm','medium','stress']:
        x=net.query('asset_i==@i and asset_j==@j and regime==@s').iloc[0]
        rr.append(f'{x.n_both}/{x.n_i} ({100*x.finite_q_cond_j_given_i:.1f}\\%)')
    rows.append(rr)
tab('network','lccc',[r'Condition $\to$ target & Calm & Medium & Stress \\'],rows)
rows=[]
for d in [6,7]:
    data=rd(U[d],'covar_mc_summary.csv').query("date=='2024-01-02'")
    for r in data.itertuples():
        rows.append([d,r.condition+' $\\to$ '+r.target,fx(100*r.CoVaR,3),fx(100*r.CoES,3),fx(100*r.CoVaR_mc_sd,3),fx(100*r.CoES_mc_sd,3),r.min_joint_tail])
tab('covar','rlrrrrr',[r'$d$ & Condition $\to$ target & CoVaR (\%) & CoES (\%) & VaR SD (pp) & ES SD (pp) & Min. tail $n$ \\'],rows)
# Prespecified model-pair uncertainty intervals under varying block length.
rows=[]
for d in [6,7]:
    p=rd(U[d],'paired_inference.csv')
    for m in ['rolling_matched','WRR_no_tail','WRR_resampled_B20']:
        a=p.query("benchmark==@m and metric=='nll'").set_index('block')
        rows.append([d,LAB[m],fx(a.difference.iloc[0],4),*[f'[{a.loc[l,"ci_low"]:.4f}, {a.loc[l,"ci_high"]:.4f}]' for l in [10,20,40]],pfx(a.p_holm.iloc[0])])
tab('blocks','rlrcccr',[r'$d$ & Comparator & Diff. & Block 10 & Block 20 & Block 40 & Holm $p$ \\'],rows)
# Figures: legible vector PDFs, print-safe palette.
plt.rcParams.update({'font.size':10,'axes.labelsize':10,'axes.titlesize':10,'legend.fontsize':8.5,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False,'figure.dpi':140,'savefig.bbox':'tight'})
navy='#145374'; teal='#19827b'; gold='#ad6b26'; red='#a33e40'; grey='#777777'
reg['date']=pd.to_datetime(reg.date)
fig,ax=plt.subplots(figsize=(7.5,2.65))
ax.plot(reg.date,reg.lagged_vix,c=navy,lw=.85)
cuts=rd(U[6],'regime_cutoffs.csv').cut.to_numpy()
ax.axhspan(cuts[1],reg.lagged_vix.max()+3,color=red,alpha=.07)
ax.axhline(cuts[0],ls='--',lw=1,color=teal,label=f'60th percentile = {cuts[0]:.2f}')
ax.axhline(cuts[1],ls='--',lw=1,color=red,label=f'85th percentile = {cuts[1]:.2f}')
for date in ['2022-01-01','2024-01-01']:ax.axvline(pd.Timestamp(date),c=grey,lw=.85,ls=':')
ax.text(pd.Timestamp('2017-06-01'),ax.get_ylim()[1]*.95,'Training',ha='center',va='top',fontsize=9)
ax.text(pd.Timestamp('2023-01-01'),ax.get_ylim()[1]*.95,'Validation',ha='center',va='top',fontsize=9)
ax.text(pd.Timestamp('2025-01-01'),ax.get_ylim()[1]*.95,'Test',ha='center',va='top',fontsize=9)
ax.set(xlabel='Date',ylabel='Lagged VIX level');ax.xaxis.set_major_locator(mdates.YearLocator(2));ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y'));ax.legend(loc='upper left',frameon=False,ncol=2,bbox_to_anchor=(0,-.23));ax.grid(axis='y',alpha=.18)
fig.savefig(F/'fig01_vix_regimes.pdf');plt.close(fig)
fig,axs=plt.subplots(2,1,figsize=(7.4,4.5),sharex=True,sharey=True)
for (d,b),ax in zip(U.items(),axs):
    x=rd(b,'daily_forecasts.csv')
    sel=x.query("model in ['WRR_direct','rolling_matched','WRR_no_tail']")
    for m,col,ls in [('WRR_direct',navy,'-'),('rolling_matched',grey,'--'),('WRR_no_tail',teal,':')]:
        z=sel[sel.model==m].sort_values('date');ax.plot(pd.to_datetime(z.date),100*z.probability,label=LAB[m].replace('$',''),c=col,ls=ls,lw=1.0)
    ref=sel.query("model=='WRR_direct'")
    ax.scatter(pd.to_datetime(ref.loc[ref.observed==1,'date']),np.full((ref.observed==1).sum(),.15),s=8,c=red,marker='|',label='Realized distress',zorder=5)
    ax.set(ylabel='Probability (%)',title=f'({chr(97+d-6)}) {d}-variable copula')
    ax.grid(axis='y',alpha=.18)
axs[-1].set(xlabel='Date');axs[-1].xaxis.set_major_locator(mdates.MonthLocator(interval=4));axs[-1].xaxis.set_major_formatter(mdates.DateFormatter('%b %Y')); fig.autofmt_xdate(rotation=0)
handles,labels=axs[0].get_legend_handles_labels();fig.legend(handles,labels,loc='upper center',ncol=4,frameon=False,bbox_to_anchor=(.5,1.02));fig.tight_layout(rect=(0,0,1,.96));fig.savefig(F/'fig02_event_probabilities.pdf');plt.close(fig)
# Exploratory conditional log-score contrasts; no significance markings.
rr=pd.read_csv(A/'exploratory_regime_inference.csv').query("spec=='primary' and block==20 and contrast!='stress_minus_calm'")
fig,axs=plt.subplots(1,2,figsize=(7.2,2.5),sharey=True)
for d,ax in zip([6,7],axs):
    z=rr.query('dimension==@d').set_index('contrast').loc[['calm','medium','stress']]
    est=z.estimate.to_numpy()
    ax.axvline(0,color=grey,ls='--',lw=.8)
    for y,(row,col) in enumerate(zip(z.itertuples(),[teal,grey,red])):
        ax.plot([row.ci_low,row.ci_high],[y,y],c=col,lw=1.7)
        ax.plot(row.estimate,y,'o',ms=6,c=col)
    ax.set(yticks=range(3),yticklabels=['Calm (241)','Medium (210)','Stress (51)'],title=f'{d}D',xlabel='WRR minus matched log score',xlim=(-.32,.46));ax.grid(axis='x',alpha=.16)
axs[0].invert_yaxis();fig.tight_layout();fig.savefig(F/'fig03_regime_logscore.pdf');plt.close(fig)
fig,axs=plt.subplots(2,1,figsize=(7.4,3.8),sharex=True,sharey=True)
for d,ax in zip([6,7],axs):
    z=unc.query('dimension==@d');dates=pd.to_datetime(z.date)
    ax.errorbar(dates,100*z.boot_mean,yerr=np.vstack((100*(z.boot_mean-z.boot_low),100*(z.boot_high-z.boot_mean))),ls='none',marker='o',capsize=3,c=navy,label='Local refits: mean and 95% interval')
    ax.set(ylabel='Probability (%)',title=f'{d}D at nine fixed origins')
    ax.grid(axis='y',alpha=.2)
axs[-1].set(xlabel='Date');axs[-1].xaxis.set_major_locator(mdates.MonthLocator(interval=4));axs[-1].xaxis.set_major_formatter(mdates.DateFormatter('%b %Y'));fig.tight_layout()
fig.savefig(F/'fig04_local_uncertainty.pdf');plt.close(fig)
print(f'Generated {len(list(T.glob("*.tex")))} tables and {len(list(F.glob("*.pdf")))} figures.')

# Final manuscript: preserve the archived score tables and add the full
# matched comparison plus the analytical tail-tilt and state diagnostics.
local_summary=pd.read_csv(A/'local_logscore_uncertainty.csv')
tlt=pd.read_csv(A/'tail_tilt_mass.csv')
cohit=pd.read_csv(A/'training_cohit_directed.csv')
joint=pd.read_csv(A/'joint6_vix_upper_event.csv')
marg=pd.read_csv(R/'tables/marginal_selection.csv').query('selected==True').set_index('asset')
tab('pit','llrrr',[r"Asset & Variance model & 5\% hits (\%) & 95\% interval (\%) & Holm $p$ \\"],
    [[r.asset,{'eGARCH':'EGARCH','gjrGARCH':'GJR-GARCH','sGARCH':'GARCH'}[marg.loc[r.asset,'model']],
      fx(100*r.hit_rate,2),f'[{100*r.hit_ci_low:.2f}, {100*r.hit_ci_high:.2f}]',pfx(r.holm_p_documented)] for r in cal.itertuples()])
# State counts include events so the scarcity of stress events is visible.
evt=pd.read_csv(R/'main_no_vix/tables/regime_event_counts.csv')
st_events=evt.query('observed==1').set_index('regime').Freq
tab('states','lrrrr',[r'Period & Calm & Medium & Stress & Total \\'],
    [[s.title(),*rgp, sum(rgp)] for s in ['train','validation','test']
        for rgp in [pd.crosstab(reg.split,reg.regime).loc[s,['calm','medium','stress']].tolist()]]
    +[['Test distress events',*[int(st_events[s]) for s in ['calm','medium','stress']],int(st_events.sum())]])
# Matched comparisons of all five inference losses.
rows=[]
metlabels={'nll':'NLL','energy_common':'Energy','variogram_common':'Variogram','brier':'Brier','quantile_loss':'Quantile loss (bp)'}
for metric in metlabels:
    rr=[metlabels[metric]]
    for d in [6,7]:
        x=P.query("spec=='primary' and dimension==@d and benchmark=='rolling_matched' and metric==@metric").iloc[0]
        scale=1e4 if metric=='quantile_loss' else 1
        dp=3 if scale!=1 else (4 if metric=='nll' else 6)
        rr += [fx(scale*x.difference,dp),f'[{scale*x.ci_low:.{dp}f}, {scale*x.ci_high:.{dp}f}]',pfx(x.p_holm)]
    rows.append(rr)
tab('matched','lrcrrcr',[r' & \multicolumn{3}{c}{6D} & \multicolumn{3}{c}{7D} \\',r'\cmidrule(lr){2-4}\cmidrule(lr){5-7}',r'Loss & Diff. & 95\% interval & $p_H$ & Diff. & 95\% interval & $p_H$ \\'],rows)
# All core NLL comparisons and main mechanism effects; matched NLL is
# given in the dedicated table.
rows=[]
for model,metric in [(m,'nll') for m in ['rolling_vine','dcc_gaussian','t_copula','static_vine','regime_vine','gaussian','WRR_raw','WRR_no_tail','WRR_resampled_B20']]+[(m,'variogram_common') for m in ['WRR_raw','WRR_no_tail']]:
    rr=[LAB[model],'NLL' if metric=='nll' else 'Variogram']
    for d in [6,7]:
        x=P.query("spec=='primary' and dimension==@d and benchmark==@model and metric==@metric").iloc[0]
        rr += [fx(x.difference,4),f'[{x.ci_low:.4f}, {x.ci_high:.4f}]',pfx(x.p_holm)]
    rows.append(rr)
tab('inference','llrcrrcr',[r' & & \multicolumn{3}{c}{6D} & \multicolumn{3}{c}{7D} \\',r'\cmidrule(lr){3-5}\cmidrule(lr){6-8}',r'Comparator & Loss & Diff. & 95\% interval & $p_H$ & Diff. & 95\% interval & $p_H$ \\'],rows)
# Direct estimates of the event tilt before ESS adjustment, and the actual
# final ESS and rank diagnostics. Equal means across 51 refits.
ww=[]
for d in [6,7]:
    di=rd(U[d],'refit_weight_diagnostics.csv')
    wr=di.query("model=='WRR_direct'")
    no=di.query("model=='WRR_no_tail'")
    en=di.query("model=='WRR_resampled_B20'")
    a=tlt.query('dimension==@d').mean(numeric_only=True)
    ww.append([str(int(wr.n_window.iloc[0])),fx(wr.ess.mean(),1),fx(100*(wr.ess/wr.n_window).mean(),2),
      fx(100*(wr.ess/wr.n_window).min(),2),str(int((wr.uniform_mix>0).sum())),fx(no.ess.mean(),1),
      fx(en.unique_mean.mean(),1),fx(a.event_mass_base*100,2),fx(a.event_mass_tilted*100,2),
      fx(a.mean_marginal_hit_base*100,2),fx(a.mean_marginal_hit_tilted*100,2),
      fx(wr.weighted_ks_before.mean(),5),fx(wr.weighted_ks_after.mean(),5)])
wl=['Window length','Mean effective sample size',r'Mean ESS / window (\%)',r'Minimum ESS / window (\%)',
    'Refits with uniform mixing','ESS without tail bonus','Distinct rows per ensemble component',
    r'Base event mass (\%)',r'Tilted event mass (\%)',r'Base mean marginal hit mass (\%)',
    r'Tilted mean marginal hit mass (\%)','Weighted KS distance before ranks','Weighted KS distance after ranks']
rows=[[l,ww[0][i],ww[1][i]] for i,l in enumerate(wl)]
tab('weights','lrr',[r'Diagnostic & 6D & 7D \\'],rows)
# Include log score uncertainty alongside event probability. Values are
# summaries of fresh local bootstrap fits, not uncertainty on the mean score.
rows=[]
for d in [6,7]:
    a=unc.query('dimension==@d');l=local_summary.query('dimension==@d')
    cs=rd(U[d],'covar_mc_summary.csv')
    rows.append([1791,1791,fx(100*(a.boot_high-a.boot_low).mean(),3),fx(100*a.boot_sd.mean(),3),
            fx(l.log_score_sd.mean(),3),fx(l.log_score_width.mean(),3),
            *[fx(100*a['mc_sd_'+str(n)].mean(),3) for n in [2000,5000,20000]],int(cs.min_joint_tail.min())])
ul=['Local refits attempted','Successful local refits','Mean probability interval width (pp)','Refit probability SD (pp)',
    'Refit log-score SD','Mean log-score interval width','Probability MC SD: $M=2{,}000$ (pp)',
    'Probability MC SD: $M=5{,}000$ (pp)','Probability MC SD: $M=20{,}000$ (pp)','Minimum conditional-tail count']
tab('uncertainty','lrr',[r'Diagnostic & 6D & 7D \\'],[[l,*[rows[j][i] for j in [0,1]]] for i,l in enumerate(ul)])
# Additional joint-stress event includes the VIX upper decile.
jt=joint.set_index('model')
tab('vix_joint','lrr',[r'Model & Mean forecast (\%) & Brier loss \\'],
    [[LAB[m],fx(100*jt.loc[m,'probability'],2),fx(jt.loc[m,'brier'],6)] for m in base+['WRR_no_tail','WRR_resampled_B20']])
# A compact sensitivity table reports all paired matched effects rather
# than ranking the absolute densities of different marginal transforms.
rows=[]
for spec in ['empirical_margin','train_winsor']:
    for metric in metlabels:
        rr=[rsp[spec],metlabels[metric]]
        for d in [6,7]:
            v=P.query("spec==@spec and dimension==@d and benchmark=='rolling_matched' and metric==@metric").iloc[0]
            scale=1e4 if metric=='quantile_loss' else 1
            dp=3 if scale!=1 else 6
            rr +=[fx(v.difference*scale,dp),pfx(v.p_holm)]
        rows.append(rr)
tab('sensitivity_matched','llrrrr',[r' & & \multicolumn{2}{c}{6D} & \multicolumn{2}{c}{7D} \\',r'\cmidrule(lr){3-4}\cmidrule(lr){5-6}',r'Specification & Loss & Diff. & $p_H$ & Diff. & $p_H$ \\'],rows)
# Heatmaps of observed training co-hit associations by state. Rows are
# conditions and columns are targets. Omit diagonals (self-conditioning).
assets=['XLE','ICLN','TAN','XLF','GLD','UUP']
fig,axs=plt.subplots(2,2,figsize=(6.7,5.6),layout='constrained')
mat={}
for st_ in ['calm','medium','stress']:
    mat[st_]=cohit.query('state==@st_').pivot(index='condition',columns='target',values='conditional').reindex(index=assets,columns=assets).to_numpy()
mat['stress minus calm']=mat['stress']-mat['calm']
for st_,ax in zip(mat,axs.flat):
    z=mat[st_]
    diff=st_=='stress minus calm'
    vmax=max(.3,np.nanmax(np.abs(z))) if diff else .8
    im=ax.imshow(np.ma.masked_invalid(z),cmap='RdBu_r' if diff else 'YlGnBu',vmin=-vmax if diff else 0,vmax=vmax,aspect='equal')
    for i in range(6):
        for j in range(6):
            if i!=j:
                val=z[i,j]
                ax.text(j,i,f'{val:+.2f}' if diff else f'{val:.2f}',ha='center',va='center',fontsize=8.3,color='white' if (abs(val)>.62*vmax if diff else val>.48) else '#172b3a')
    ax.set(xticks=np.arange(6),xticklabels=assets,yticks=np.arange(6),yticklabels=assets,title=st_.capitalize())
    ax.tick_params(axis='x',labelrotation=38,labelsize=8.2)
    ax.tick_params(axis='y',labelsize=8.2)
    ax.spines[['left','bottom']].set_visible(False);ax.tick_params(length=0)
    if diff:
        fig.colorbar(im,ax=ax,shrink=.83,pad=.025,label='Frequency difference')
    else:
        fig.colorbar(im,ax=ax,shrink=.83,pad=.025,label='Conditional frequency')
for ax in axs[:,0]:ax.set_ylabel('Conditioning asset')
for ax in axs[1,:]:ax.set_xlabel('Target asset')
fig.savefig(F/'fig05_training_cohits.pdf');plt.close(fig)
# Complete principal test-period PIT diagnostics. All entries are derived
# from the archived threshold hits and the audited reporting p-values.
pc_all=pd.read_csv(A/'pit_calibration_checked.csv').query("spec=='primary' and split=='test'")
qs=[.01,.05,.10,.50,.90,.95,.99]
rows=[]
for asset in ['XLE','ICLN','TAN','XLF','GLD','UUP','VIX']:
    a=pc_all.query('asset==@asset').set_index('q').loc[qs]
    assert a.n.eq(502).all()
    rr=[asset]
    for q,r in a.iterrows():
        v=fx(r.hit_rate,3)
        if pd.notna(r.holm_p_documented) and r.holm_p_documented<.05:
            v=r'$'+v+r'^{*}$'
        rr.append(v)
    z=a.iloc[0]
    rr.extend([fx(z.ks_stat,4),pfx(z.ks_p_iid_reference),pfx(z.lb_z_p),pfx(z.lb_z2_p)])
    rows.append(rr)
tab('pit_full','lrrrrrrrrrrr',
    [r' & \multicolumn{7}{c}{PIT hit rate at nominal threshold $q$} & \multicolumn{2}{c}{KS} & \multicolumn{2}{c}{LB $p$} \\',
     r'\cmidrule(lr){2-8}\cmidrule(lr){9-10}\cmidrule(lr){11-12}',
     r'Asset & 0.01 & 0.05 & 0.10 & 0.50 & 0.90 & 0.95 & 0.99 & $D_{\rm KS}$ & $p_{\rm KS}$ & $\eta$ & $\eta^2$ \\'],rows)
print('Final assets:',len(list(T.glob('*.tex'))),'tables,',len(list(F.glob('*.pdf'))),'figures')
