#!/usr/bin/env python3
"""Independent, read-only audit of the completed WRR-TVine CSV archive.

Run: python audit_results.py extracted/results/full output_dir
Requires Python 3, numpy, pandas, scipy. No model refits, downloads or tuning.
The optional regime contrasts are retrospective exploratory analysis.
"""
from pathlib import Path
import sys, json
import numpy as np
import pandas as pd
from scipy.stats import norm, kendalltau
from scipy.optimize import brentq

ROOT=Path(sys.argv[1]); OUT=Path(sys.argv[2]); OUT.mkdir(parents=True,exist_ok=True)
A=['XLE','ICLN','TAN','XLF','GLD','UUP','VIX']
U_NAMES=['main_no_vix','extended_with_vix']
SPECS=['primary','empirical_margin','train_winsor']
METRICS=['nll','energy_common','variogram_common','brier','quantile_loss']
MODELS={'gaussian','t_copula','static_vine','regime_vine','dcc_gaussian','WRR_direct','rolling_vine','rolling_matched'}
EXTRA={'WRR_raw','WRR_no_tail','WRR_no_recency','WRR_no_regime','WRR_no_vix_distance','WRR_other_family','WRR_legacy_resampled_B1','WRR_resampled_B1','WRR_resampled_B5','WRR_resampled_B20','WRR_threshold_55_80','WRR_threshold_65_90'}
CHECKS=[]
def ck(ok,label):
    if not bool(ok):
        CHECKS.append(dict(check=label,result='FAIL'))
        print('FAIL',label,flush=True)
    else:
        CHECKS.append(dict(check=label,result='PASS'))
def eq(x,y,label,atol=5e-10,rtol=1e-9):
    ck(np.allclose(x,y,atol=atol,rtol=rtol,equal_nan=True),label)
def rd(base,name):
    return pd.read_csv(base/'tables'/name,low_memory=False)
def specpath(spec):
    return ROOT if spec=='primary' else ROOT/'sensitivity'/spec
def hac(d,lag=10):
    d=np.asarray(d,float); n=len(d); z=d-d.mean()
    v=np.dot(z,z)/n
    for j in range(1,min(lag,n-2)+1): v+=2*(1-j/(lag+1))*np.dot(z[j:],z[:-j])/n
    se=np.sqrt(max(0,v)/n)
    p=2*norm.sf(abs(d.mean()/se)) if se>1e-15 else (1. if abs(d.mean())<1e-15 else 0.)
    return se,p
def holm(p):
    p=np.asarray(p,float); order=np.argsort(p,kind='stable'); adj=np.minimum(1,np.maximum.accumulate(p[order]*np.arange(len(p),0,-1))); out=np.empty(len(p)); out[order]=adj; return out
def block_indices(n,l,B,seed):
    r=np.random.default_rng(seed)
    return ((r.integers(0,n,size=(B,(n+l-1)//l))[:,:,None]+np.arange(l))%n).reshape(B,-1)[:,:n]
def tied_mid_ks(X,w):
    b=[]; a=[]
    for col in X.T:
        v,g=np.unique(col,return_inverse=True)
        mass=np.bincount(g,w); c=np.cumsum(mass); mids=np.clip(c-mass/2,1e-8,1-1e-8)
        b.append(np.max(np.maximum(abs(c-v),abs(c-mass-v))))
        a.append(np.max(np.maximum(abs(c-mids),abs(c-mass-mids))))
    return np.mean(b),np.mean(a)
prices=rd(ROOT,'prices_snapshot.csv').set_index('date')
raw=rd(ROOT,'untrimmed_returns.csv').set_index('date')
pit=rd(ROOT,'marginal_PIT.csv').set_index('date')
dreturns=np.log(prices[A]).diff().iloc[1:]
ck(len(prices)==4024 and prices.index.is_unique and prices.index.is_monotonic_increasing,'4024 sorted unique common prices')
ck(np.isfinite(prices[A]).all().all() and (prices[A]>0).all().all(),'finite positive common prices')
ck(len(dreturns)==4023 and len(raw)==len(pit)==3273,'return counts and 750-day burn-in')
ck(pit.index.equals(raw.index) and pit.index.equals(dreturns.index[750:]),'date alignment after burn-in')
eq(raw[A],dreturns.loc[raw.index,A],'untrimmed log returns reproduce frozen prices')
ck(pit.split.value_counts().to_dict()=={'train':2270,'test':502,'validation':501},'chronological sample sizes')
ck(pit.index[pit.split=='train'][-1]=='2021-12-31' and pit.index[pit.split=='test'][0]=='2024-01-02','split boundaries')
ck((pit[A]>0).all().all() and (pit[A]<1).all().all(),'PITs strictly interior')
ck(len(rd(ROOT,'dropped_prices.csv'))==0,'no invalid prices in aligned intersection')
test=pit.index[pit.split=='test']; train=pit.split=='train'
obsportfolio=np.expm1(raw.loc[test,A[:6]]).mean(axis=1)
v=prices.VIX.shift(1).loc[pit.index].to_numpy()
cut=np.quantile(v[train],[.60,.85],method='median_unbiased')
states=np.where(v<=cut[0],'calm',np.where(v<=cut[1],'medium','stress'))
ep=pd.DataFrame(index=pit.index,columns=A,dtype=float)
R=dreturns.to_numpy()
for ix,t in enumerate(range(750,len(R))):
    past=R[t-750:t]; ep.iloc[ix]=((past<R[t]).sum(axis=0)+.5*(past==R[t]).sum(axis=0)+.5)/751
PITS={'primary':pit[A],'empirical_margin':ep}
selected=[]; pit_corrections=[]; cases={}; regimes=[]; grid_summary=[]; scores=[]; extra_stress=[]; refit_details=[]; uncertainty=[]; comparisons=[]
for spec in SPECS:
    base=specpath(spec)
    calib=rd(base,'PIT_calibration.csv')
    ck(len(calib)==98 and not calib.duplicated(['asset','split','q']).any(),f'{spec}: 98 unique PIT diagnostics')
    after=calib.copy()
    deg=(after.hit_rate==0)|(after.hit_rate==1)
    after['degenerate_hit_series']=deg
    after['hac_p_documented']=after.calibration_hac_p.mask(deg)
    after['holm_p_documented']=np.nan
    after['spec']=spec
    for s,group in after.groupby('split'):
        p=group.hac_p_documented.fillna(1)
        after.loc[group.index,'holm_p_documented']=pd.Series(holm(p),index=group.index).mask(group.degenerate_hit_series)
    pit_corrections.append(after)
    if spec in PITS:
        for row in calib.itertuples():
            ui=PITS[spec].loc[pit.split==row.split,row.asset]
            hit=(ui.to_numpy()<=row.q).astype(float)
            eq(hit.mean(),row.hit_rate,f'{spec}: PIT {row.asset}/{row.split}/{row.q:g}')
    if spec!='empirical_margin':
        m=rd(base,'marginal_selection.csv')
        ck(len(m)==21 and m.selected.sum()==7 and (m.final_convergence==0).all() and np.isfinite(m.validation_nll).all(),f'{spec}: 21 usable GARCH candidates')
        for a,g in m.groupby('asset'):
            ck(g.loc[g.selected,'model'].iloc[0]==g.loc[g.validation_nll.idxmin(),'model'],f'{spec}: margin validation minimum {a}')
    for un in U_NAMES:
        d=6 if un=='main_no_vix' else 7; name=f'{spec}/{d}D'; b=base/un
        reg=rd(b,'regimes.csv').set_index('date'); rc=rd(b,'regime_cutoffs.csv')
        ck(reg.index.equals(pit.index),name+': state calendar')
        eq(reg.lagged_vix,v,name+': lagged VIX reproduces prices')
        eq(rc.cut,cut,name+': training-only regime cutoffs')
        ck(np.array_equal(reg.regime,states),name+': regime labels')
        x=rd(b,'daily_forecasts.csv'); g=rd(b,'tail_grid_daily.csv')
        want=MODELS|EXTRA if spec=='primary' else MODELS
        ck(set(x.model)==want and len(x)==502*len(want),name+': expected daily rows and model set')
        ck(not x.duplicated(['model','date']).any() and np.isfinite(x.select_dtypes('number')).all().all(),name+': no missing/duplicate daily forecasts')
        ck(((x.probability>=0)&(x.probability<=1)&x.observed.isin([0,1])).all(),name+': valid probabilities and hits')
        ck(x.groupby('model').date.apply(lambda y: sorted(y)==list(test)).all(),name+': all models share the complete test calendar')
        ck((x.regime==x.date.map(reg.regime)).all(),name+': forecast state alignment')
        eq(x.portfolio_return,x.date.map(obsportfolio),name+': realized six-ETF portfolio')
        eq(x.brier,(x.probability-x.observed)**2,name+': daily Brier')
        eq(x.brier_mc_adjusted,x.brier-x.probability*(1-x.probability)/9999,name+': MC Brier correction')
        eq(x.probability_mc_se,np.sqrt(x.probability*(1-x.probability)/10000),name+': MC probability SE')
        eq(x.censored_nll,np.where(x.observed==1,x.nll,-np.log1p(-np.clip(x.probability,1e-10,1-1e-10))),name+': censored log loss')
        eq(x.violation,(x.portfolio_return<=x.VaR).astype(int),name+': VaR violations')
        eq(x.quantile_loss,(.05-(x.portfolio_return<x.VaR).astype(int))*(x.portfolio_return-x.VaR),name+': quantile loss')
        eq(x.es_identification,x.ES-x.VaR+(x.VaR-x.portfolio_return)*x.violation/.05,name+': ES identification')
        ck((x.ES<=x.VaR+1e-14).all() and (x.n_es_draws>=500).all(),name+': tail risk order and simulation support')
        ck((x.groupby('date')['observed'].nunique()==1).all(),name+': common event across models')
        ck(len(g)==6*len(x) and not g.duplicated(['date','model','q','k']).any() and len(g[['q','k']].drop_duplicates())==6,name+': all six tail event grids')
        eq(g.brier,(g.probability-g.observed)**2,name+': tail-grid daily Brier')
        ck(((g.probability>=0)&(g.probability<=1)&g.observed.isin([0,1])).all(),name+': tail grid bounds')
        if spec in PITS:
            U=PITS[spec]; ev=((U.loc[test,A[:6]]<=.1).sum(axis=1)>=3).astype(int)
            eq(x.observed,x.date.map(ev),name+': main observed event reproduces PITs')
            for (q,k),group in g.groupby(['q','k']):
                e=(U.loc[test,A[:6]].le(q).sum(axis=1)>=k).astype(int)
                eq(group.observed,group.date.map(e),name+f': tail event {q:g}/{k}')
            if d==7:
                joint=ev*(U.loc[test,'VIX']>=.9).astype(int)
                eq(x.joint6_vix_upper_observed,x.date.map(joint),name+': observed VIX-up joint event')
                ck((x.joint6_vix_upper_probability<=x.probability+1e-13).all(),name+': VIX-up subset probability')
        # Model-average and conditional summaries.
        s=rd(b,'model_comparison.csv').set_index('model').sort_index()
        means=x.groupby('model')[s.columns.drop(['n','mean_log_score']).tolist()].mean().sort_index()
        eq(s[means.columns],means,name+': all stored score means')
        eq(s.mean_log_score,-s.nll,name+': log score sign')
        eq(s.n,x.groupby('model').size().sort_index(),name+': summary n')
        br=rd(b,'by_regime.csv').set_index(['model','regime']).sort_index()
        eq(br,x.groupby(['model','regime'])[br.columns].mean().sort_index(),name+': by-regime summary means')
        gg=g.groupby(['model','q','k'],as_index=False)[['brier','probability','observed']].mean()
        gg.insert(0,'dimension',d); gg.insert(0,'spec',spec); grid_summary.append(gg)
        ss=s.reset_index();ss.insert(0,'dimension',d);ss.insert(0,'spec',spec);scores.append(ss)
        # Independent forecasts loss-differences, Bartlett HAC and Holm.
        f=rd(b,'paired_inference.csv')
        ck(len(f)==(len(want)-1)*5*3 and set(f.block)=={10,20,40} and (f.bootstrap_B==1999).all() and (f.n==502).all(), name+': inference families, 1999 draws, block sizes')
        w=x.set_index(['model','date']).sort_index()
        ref=w.loc['WRR_direct'].loc[test]
        for (m,metric),group in f.groupby(['benchmark','metric']):
            diff=(ref[metric]-w.loc[m].loc[test,metric]).to_numpy()
            se,p=hac(diff)
            eq(group.difference,diff.mean(),name+f': paired mean {m}/{metric}')
            eq(group.hac_se,se,name+f': HAC SE {m}/{metric}')
            eq(group.p_two_sided,p,name+f': HAC p {m}/{metric}')
        for l,group in f.groupby('block'):
            eq(group.p_holm,holm(group.p_two_sided),name+f': Holm {l}')
        ck((f.ci_low<=f.ci_high).all(),name+': ordered bootstrap intervals')
        f20=f[f.block==20].copy();f20.insert(0,'dimension',d);f20.insert(0,'spec',spec);comparisons.append(f20)
        # Validation-selected hyperparameters and rolling refit calendars.
        di=rd(b,'refit_weight_diagnostics.csv')
        sel=rd(b,'selected_wrr.csv').iloc[0]
        for var in ['wrr','rolling']:
            t=rd(b,f'tuning_{var}.csv')
            best=t.sort_values(['validation_nll','window','refit','family']).iloc[0]
            ck(len(t)==12 and (t.n_validation==501).all() and np.isfinite(t.validation_nll).all(),name+f': 12 complete {var} validation candidates')
            ts=sel if var=='wrr' else best
            ck(all(ts[k]==best[k] for k in ['window','refit','family']),name+f': {var} grid minimum')
            mdl='WRR_direct' if var=='wrr' else 'rolling_vine'
            dr=di[di.model==mdl]
            ck(list(dr.t)==list(range(2772,3274,int(best.refit))) and (dr.n_window==best.window).all(),name+f': {var} refit schedule')
            selected.append(dict(spec=spec,dimension=d,model=mdl,window=best.window,refit=best.refit,family=best.family,validation_log_score=-best.validation_nll))
        ck(di.model.isin(want).all() and np.isfinite(di.drop(columns=['unique_mean']).select_dtypes('number')).all().all() and (di.ess>=.15*di.n_window-.02).all(),name+': ESS support and refit diagnostics')
        dd=rd(b,'dcc_validation.csv')
        ck(len(dd)==12 and (dd.a+dd.b<.995).all() and np.isfinite(dd.validation_nll).all(),name+': DCC grid')
        # Rebuild Gaussian and DCC predictive log densities directly from archived PITs.
        if spec in PITS:
            U=PITS[spec].iloc[:,:d].to_numpy()
            z=norm.ppf(U); S=np.corrcoef(z[train],rowvar=False)
            def nll(z,R):
                return .5*(np.linalg.slogdet(R)[1]+np.einsum('ij,jk,ik->i',z,np.linalg.inv(R)-np.eye(d),z))
            eq(w.loc['gaussian'].loc[test,'nll'],nll(z[pit.split=='test'],S),name+': Gaussian predictive density rebuilt',atol=3e-8)
            losses=[]
            for t in dd.itertuples():
                Q=S.copy(); lls=[]
                for i in range(len(z)):
                    if i: Q=(1-t.a-t.b)*S+t.a*np.outer(z[i-1],z[i-1])+t.b*Q
                    if i>=2270:
                        sd=np.sqrt(np.diag(Q)); cor=Q/np.outer(sd,sd); lls.append(nll(z[i:i+1],cor)[0])
                losses.append(lls)
            losses=np.array(losses)
            eq(dd.validation_nll,losses[:,:501].mean(axis=1),name+': all DCC validation losses rebuilt',atol=3e-8)
            eq(w.loc['dcc_gaussian'].loc[test,'nll'],losses[np.argmin(dd.validation_nll),501:],name+': DCC predictive density rebuilt',atol=3e-8)
            # Every recorded rolling weight is reconstructed; bootstrap prefix
            # uniqueness is deliberately not asserted because those draws are not archived.
            errors=[]; per=[]
            for t in di.itertuples():
                start=t.t-1-int(t.n_window); stop=t.t-1; past=slice(start,stop)
                arr=U[past]; vs=v[past]; n=len(arr)
                state=states
                if t.model.startswith('WRR_threshold_'):
                    pp=[int(i)/100 for i in t.model.split('_')[-2:]]
                    ct=np.quantile(v[train],pp,method='median_unbiased')
                    state=np.where(v<=ct[0],'calm',np.where(v<=ct[1],'medium','stress'))
                if t.model in ['rolling_vine','rolling_matched']:
                    p=np.ones(n)/n; rawess=float(n); mix=0.
                else:
                    ll=-(0 if t.model=='WRR_no_recency' else 2)*np.arange(n-1,-1,-1)/n
                    ll+=np.log(np.where(state[past]==state[stop],1 if t.model=='WRR_no_regime' else 3,1))
                    ll-=(0 if t.model=='WRR_no_vix_distance' else 1.25)*abs(vs-v[stop])/np.std(vs,ddof=1)
                    ll+=np.log1p((0 if t.model=='WRR_no_tail' else 2)*((arr[:,:6]<=.1).sum(axis=1)>=3))
                    p=np.exp(ll-ll.max());p/=p.sum();rawess=1/(p@p)
                    mix=t.uniform_mix # use recorded uniroot mixing coefficient, check floor separately
                    p=(1-mix)*p+mix/n
                before,afterks=tied_mid_ks(arr,p)
                if t.model in ['WRR_raw','WRR_legacy_resampled_B1','rolling_vine']: afterks=before
                errors.append([abs(rawess-t.raw_ess),abs(1/(p@p)-t.ess),abs(max(p)-t.max_weight),abs(before-t.weighted_ks_before),abs(afterks-t.weighted_ks_after)])
                if t.model=='WRR_direct': per.append(dict(date=t.date,ess_fraction=t.ess/t.n_window,ks_before=before,ks_after=afterks,regime=state[stop]))
            errs=np.array(errors).max(axis=0)
            ck((errs<1e-7).all(),name+': reconstructed all rolling weights and weighted KS')
            if (errs>=1e-7).any(): print('max weight reconstruction errors',errs)
            for r in per: r.update(spec=spec,dimension=d);refit_details.append(r)
            net=rd(b,'finite_threshold_network.csv')
            vals=[]
            for r in net.itertuples():
                vv=PITS[spec].loc[train & (states==r.regime),[r.asset_i,r.asset_j]].to_numpy(); hit=vv<=.1
                vals.append([len(vv),hit[:,0].sum(),hit[:,1].sum(),np.all(hit,axis=1).sum(),kendalltau(vv[:,0],vv[:,1]).statistic])
            eq(net[['n','n_i','n_j','n_both','kendall_tau']],vals,name+': observed training network counts and Kendall tau')
        # Exploratory, calendar-preserving regime contrasts (held-out period already observed).
        delta=(w.loc['rolling_matched'].loc[test,'nll']-ref.nll).to_numpy()
        st=ref.regime.to_numpy()
        means={s:delta[st==s].mean() for s in ['calm','medium','stress']}
        nst={s:(st==s).sum() for s in means}
        g0={s:(st==s)*(delta-means[s])/(nst[s]/len(st)) for s in means}
        for contrast in ['calm','medium','stress','stress_minus_calm']:
            isdiff=contrast=='stress_minus_calm'
            value=means['stress']-means['calm'] if isdiff else means[contrast]
            influence=g0['stress']-g0['calm'] if isdiff else g0[contrast]
            se,_=hac(influence); p=2*norm.sf(abs(value)/se) if se>0 else np.nan
            for l in [10,20,40]:
                ix=block_indices(len(st),l,1999,20261125+100*d+7*SPECS.index(spec)+l)
                with np.errstate(divide='ignore',invalid='ignore'):
                    def avg(s):
                        indic=st[ix]==s
                        return (delta[ix]*indic).sum(axis=1)/indic.sum(axis=1)
                    sims=avg('stress')-avg('calm') if isdiff else avg(contrast)
                valid=np.isfinite(sims)
                cl,ch=np.quantile(sims[valid],[.025,.975])
                extra_stress.append(dict(spec=spec,dimension=d,contrast=contrast,estimate=value,n=nst['stress'] if isdiff else nst[contrast],calm_n=nst['calm'],hac_se=se,p=p,block=l,valid_draws=int(valid.sum()),ci_low=cl,ci_high=ch))
        for s,mean in means.items():regimes.append(dict(spec=spec,dimension=d,regime=s,n=nst[s],events=int(ref.loc[ref.regime==s,'observed'].sum()),WRR_minus_matched_log=mean))
        cases[spec,d]=w
        if spec=='primary':
            # Nine fixed-date origin checks, all local estimation draws, and CoVaR/CoES summaries.
            draws=rd(b,'local_estimation_bootstrap_draws.csv');ints=rd(b,'local_estimation_intervals.csv')
            mc=rd(b,'scenario_mc_replicates.csv'); co=rd(b,'covar_mc_replicates.csv'); cs=rd(b,'covar_mc_summary.csv')
            origins=test[np.round(np.linspace(1,502,9)).astype(int)-1]
            ck(len(draws)==1791 and len(ints)==9 and len(mc)==270 and len(co)==180 and len(cs)==36,name+': all uncertainty replications received')
            ck(set(draws.date)==set(origins) and set(mc.date)==set(origins) and set(co.date)==set(origins),name+': nine equally spaced origins')
            ck(not draws.duplicated(['date','b']).any() and draws.groupby('date').b.apply(lambda z:set(z)==set(range(1,200))).all() and draws[['probability','log_score']].notna().all().all(),name+': complete successful 199-block refits per date')
            for row in ints.itertuples():
                p=draws.loc[draws.date==row.date,'probability']; m=mc[mc.date==row.date]
                eq(np.quantile(p,[.025,.975]),[row.probability_low,row.probability_high],name+f': local uncertainty limits {row.date}')
                ck(row.B==row.success==199 and row.failure_rate==0 and row.interval_reportable,name+f': local success {row.date}')
                uncertainty.append(dict(dimension=d,date=row.date,regime=reg.loc[row.date,'regime'],boot_mean=p.mean(),boot_sd=p.std(),boot_low=row.probability_low,boot_high=row.probability_high,**{f'mc_sd_{M}':m.loc[m.M==M,'probability'].std() for M in [2000,5000,20000]}))
            ck(mc.groupby(['date','M']).rep.nunique().eq(10).all() and set(mc.M)=={2000,5000,20000},name+': all MC simulation repeats')
            eq(mc.brier,(mc.probability-mc.observed)**2,name+': MC diagnostic Brier')
            ck(co.reliable.all() and (co.n_joint_tail>=100).all() and co.groupby(['date','condition','target']).rep.nunique().eq(5).all(),name+': CoVaR/CoES tail support')
            for row in cs.itertuples():
                vv=co[(co.date==row.date)&(co.condition==row.condition)&(co.target==row.target)]
                eq([vv.CoVaR.mean(),vv.CoES.mean(),vv.CoVaR.std(),vv.CoES.std(),vv.n_joint_tail.min(),len(vv)],[row.CoVaR,row.CoES,row.CoVaR_mc_sd,row.CoES_mc_sd,row.min_joint_tail,row.repeats],name+f': conditional risk means and Monte Carlo SD {row.date}/{row.condition}/{row.target}')
# Cross-universe same-target inference.
cross=rd(ROOT,'common_six_cross_universe.csv')
a=cases['primary',6];b=cases['primary',7]
ck(len(cross)==80 and (cross.n==502).all(), 'cross-universe: 80 matched hypotheses')
for r in cross.itertuples():
    diff=(b.loc[r.model].loc[test,r.metric]-a.loc[r.model].loc[test,r.metric]).to_numpy();p=hac(diff)[1]
    eq([r.extended_minus_main,r.p],[diff.mean(),p],f'cross-universe: {r.model}/{r.metric}')
eq(cross.p_holm,holm(cross.p),'cross-universe Holm family')
# Observations are the same within each marginal specification; never compare
# NLL levels across dimensions or different marginal transformations.
for spec in SPECS:
    eq(cases[spec,6].loc['WRR_direct'].loc[test,'observed'],cases[spec,7].loc['WRR_direct'].loc[test,'observed'],spec+': same six-ETF event across dimensions')
stress=pd.DataFrame(extra_stress)
for l,grp in stress.groupby('block'):stress.loc[grp.index,'holm_p']=holm(grp.p.fillna(1))
report=pd.DataFrame(CHECKS)
report.to_csv(OUT/'audit_checks.csv',index=False)
pd.concat(pit_corrections,ignore_index=True).to_csv(OUT/'pit_calibration_checked.csv',index=False)
pd.DataFrame(selected).to_csv(OUT/'validation_selection.csv',index=False)
pd.concat(scores,ignore_index=True).to_csv(OUT/'score_summary.csv',index=False)
pd.concat(grid_summary,ignore_index=True).to_csv(OUT/'tail_grid_summary.csv',index=False)
pd.concat(comparisons,ignore_index=True).to_csv(OUT/'paired_block20_summary.csv',index=False)
pd.DataFrame(regimes).to_csv(OUT/'regime_logscore_contrasts.csv',index=False)
stress.to_csv(OUT/'exploratory_regime_inference.csv',index=False)
pd.DataFrame(refit_details).to_csv(OUT/'reconstructed_wrr_weights.csv',index=False)
pd.DataFrame(uncertainty).to_csv(OUT/'uncertainty_summary.csv',index=False)
print('CHECKS',report.result.value_counts().to_dict())
print('DEGENERATE HIT SERIES')
print(pd.concat(pit_corrections).query('degenerate_hit_series')[['spec','asset','split','q','hit_rate']].to_string(index=False))
print('EXPLORATORY REGIME CONTRASTS, block 20')
print(stress.query('block==20')[['spec','dimension','contrast','estimate','ci_low','ci_high','n','p','holm_p','valid_draws']].to_string(index=False))
if (report.result=='FAIL').any(): sys.exit(1)
