#!/usr/bin/env python3
"""Additional reporting from frozen daily forecasts and PITs (no model refitting).
Run: python extended_diagnostics.py <results> <analysis_dir>
Every diagnostic is conditional on the supplied estimation design.
"""
from pathlib import Path
import sys
import numpy as np
import pandas as pd

root,out=Path(sys.argv[1]),Path(sys.argv[2])
out.mkdir(parents=True,exist_ok=True)
assets=['XLE','ICLN','TAN','XLF','GLD','UUP']
pits=pd.read_csv(root/'tables/marginal_PIT.csv')
train=pits.split=='train'
counts=(pits[assets]<=.10).sum(axis=1).to_numpy(float)
events=(counts>=3).astype(float)
rows,conditional,local,refresh,cohit,tail=[] ,[] ,[] ,[] ,[] ,[]
for dim,universe in [(6,'main_no_vix'),(7,'extended_with_vix')]:
    r=root/universe/'tables'
    st=pd.read_csv(r/'regimes.csv')
    assert np.array_equal(st.date,pits.date)
    di=pd.read_csv(r/'refit_weight_diagnostics.csv')
    dr=di.query("model=='WRR_direct'").set_index('date')
    no=di.query("model=='WRR_no_tail'").set_index('date')
    v=st.lagged_vix.to_numpy(float);s=st.regime.to_numpy()
    for date,row in dr.iterrows():
        end=int(row.t)-1
        h=np.arange(end-int(row.n_window),end)
        W=len(h)
        # The original three non-tail weighting factors, before the ESS blend.
        lp=-2*(end-1-h)/W+np.log(np.where(s[h]==s[end],3,1))-1.25*abs(v[h]-v[end])/np.std(v[h],ddof=1)
        b=np.exp(lp-lp.max());b/=b.sum()
        p=b* (1+2* events[h]);p/=p.sum()
        b_final=(1-no.loc[date,'uniform_mix'])*b+no.loc[date,'uniform_mix']/W
        p_final=(1-row.uniform_mix)*p+row.uniform_mix/W
        rate0=b @ events[h];rate1=p @ events[h]
        hits0=b @ counts[h]/6;hits1=p @ counts[h]/6
        formula_rate=3*rate0/(1+2*rate0)
        formula_hits=2*((b @ (counts[h]*events[h]))-(b @ counts[h])*rate0)/(1+2*rate0)/6
        assert abs(rate1-formula_rate)<1e-12
        assert abs(hits1-hits0-formula_hits)<1e-12
        assert abs(row.ess-1/(p_final@p_final))<1e-6
        assert abs(no.loc[date,'ess']-1/(b_final@b_final))<1e-6
        rows.append(dict(dimension=dim,date=date,regime=s[end],window=W,
            event_mass_base=rate0,event_mass_tilted=rate1,mean_marginal_hit_base=hits0,
            mean_marginal_hit_tilted=hits1,mean_marginal_hit_increase=hits1-hits0,
            event_mass_no_tail_after_ess=b_final@events[h],event_mass_wrr_after_ess=p_final@events[h],
            mean_marginal_hit_no_tail_after_ess=b_final@counts[h]/6,
            mean_marginal_hit_wrr_after_ess=p_final@counts[h]/6,
            wrr_ess=row.ess,no_tail_ess=no.loc[date,'ess'],wrr_uniform_mix=row.uniform_mix,
            no_tail_uniform_mix=no.loc[date,'uniform_mix']))
    d=pd.read_csv(r/'daily_forecasts.csv')
    ref=d.query("model=='WRR_direct'").copy().sort_values('date')
    match=d.query("model=='rolling_matched'").copy().set_index('date').loc[ref.date]
    ref['matched_nll']=match.nll.to_numpy()
    ref['wrr_minus_matched_log']=ref.matched_nll-ref.nll
    rf=pd.DataFrame({'refit_date':dr.index.astype(str)})
    rf['date']=pd.to_datetime(rf.refit_date)
    ref['day']=pd.to_datetime(ref.date)
    # Forecast fits carried forward for ten common trading dates.
    rfs=pd.merge_asof(ref[['date','day','regime','wrr_minus_matched_log','observed']],rf.rename(columns={'date':'day'}),on='day',direction='backward')
    rfs['refit_regime']=rfs.refit_date.map(st.set_index('date').regime)
    assert rfs.refit_regime.notna().all()
    for (state,last),g in rfs.groupby(['regime','refit_regime']):
        refresh.append(dict(dimension=dim,state=state,refit_state=last,n=len(g),
                            mean_logscore_difference=g.wrr_minus_matched_log.mean(),events=g.observed.sum()))
    for model in ['WRR_direct','rolling_matched','WRR_no_tail','WRR_resampled_B20','dcc_gaussian']:
        md=d[d.model==model]
        for state,g in md.groupby('regime'):
            conditional.append(dict(dimension=dim,model=model,state=state,n=len(g),events=int(g.observed.sum()),
                observed=g.observed.mean(),probability=g.probability.mean(),brier=g.brier.mean(),
                adjusted_brier=g.brier_mc_adjusted.mean(),mean_log=-g.nll.mean()))
    b=pd.read_csv(r/'local_estimation_bootstrap_draws.csv')
    for date,g in b.groupby('date'):
        assert len(g)==199 and g[['log_score','probability']].notna().all().all()
        l,h=g.log_score.quantile([.025,.975])
        pl,ph=g.probability.quantile([.025,.975])
        local.append(dict(dimension=dim,date=date,regime=st.set_index('date').loc[date,'regime'],
            log_score_mean=g.log_score.mean(),log_score_sd=g.log_score.std(),log_score_low=l,log_score_high=h,
            log_score_width=h-l,probability_mean=g.probability.mean(),probability_sd=g.probability.std(),
            probability_low=pl,probability_high=ph))
    if dim==7:
        for model,g in d.groupby('model'):
            tail.append(dict(model=model,observed=g.joint6_vix_upper_observed.mean(),
                probability=g.joint6_vix_upper_probability.mean(),
                brier=np.mean((g.joint6_vix_upper_observed-g.joint6_vix_upper_probability)**2)))
    if dim==6:
        net=pd.read_csv(r/'finite_threshold_network.csv')
        for state,g in net.groupby('regime'):
            for row in g.itertuples():
                cohit.append(dict(state=state,condition=row.asset_i,target=row.asset_j,n_condition=row.n_i,n_joint=row.n_both,conditional=row.finite_q_cond_j_given_i))
                cohit.append(dict(state=state,condition=row.asset_j,target=row.asset_i,n_condition=row.n_j,n_joint=row.n_both,conditional=row.finite_q_cond_i_given_j))
        # Verify selected condition counts directly from training PITs (both directions).
        for row in cohit:
            m=train & (st.regime==row['state'])
            a=pits.loc[m,row['condition']]<=.1;c=pits.loc[m,row['target']]<=.1
            assert int(a.sum())==row['n_condition']
            assert int((a&c).sum())==row['n_joint']
            assert np.isclose((a&c).sum()/a.sum(),row['conditional'])
for name,data in [('tail_tilt_mass',rows),('state_event_forecasts',conditional),('local_logscore_uncertainty',local),
                  ('forecast_refit_state_alignment',refresh),('training_cohit_directed',cohit),('joint6_vix_upper_event',tail)]:
    pd.DataFrame(data).to_csv(out/(name+'.csv'),index=False,float_format='%.15g')
print('TAIL WEIGHT MASS (pre ESS blend)')
t=pd.DataFrame(rows)
print(t.groupby('dimension')[['event_mass_base','event_mass_tilted','mean_marginal_hit_base','mean_marginal_hit_tilted',
     'event_mass_no_tail_after_ess','event_mass_wrr_after_ess','mean_marginal_hit_no_tail_after_ess','mean_marginal_hit_wrr_after_ess']].mean().round(6).to_string())
print('STRESS EVENT CALIBRATION')
print(pd.DataFrame(conditional).query("state=='stress'").round(6).to_string(index=False))
print('LOCAL LOG SCORE UNCERTAINTY')
print(pd.DataFrame(local).groupby('dimension')[['log_score_sd','log_score_width']].mean())
print('REFIT/REALIZED STATE AGREEMENT')
print(pd.DataFrame(refresh).query("dimension==6").to_string(index=False))
print('JOINT UPPER-VIX EVENT')
print(pd.DataFrame(tail).query("model in ['WRR_direct','rolling_matched','dcc_gaussian','rolling_vine','WRR_no_tail']").round(6).to_string(index=False))
