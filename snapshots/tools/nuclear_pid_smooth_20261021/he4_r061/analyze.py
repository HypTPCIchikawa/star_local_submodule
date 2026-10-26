#!/usr/bin/env python3
import collections,csv,math,pathlib,shutil,statistics,sys
HERE=pathlib.Path(__file__).resolve().parent;sys.path.insert(0,str(HERE));import run
sys.path.insert(0,str(HERE.parent));import summarize as shared
R,O=run.ROOT,run.OUT;old=shared.load_old();cfg=run.cfg
metrics=[];atlas=[];details=[];aggregates=[];selected={};settings={};audits={};native=[]
def finite(x):return isinstance(x,(float,int)) and math.isfinite(x)
def valid(r):return r['numericalValid'] and finite(r['dataDeviance']) and finite(r['pearson']) and r['conditionalNdf']>0
raw={}
for r in run.table(O/'candidate_gate_counts/gate_counts.tsv'):raw[int(r['rBin']),r['species']]={k:float(r[k]) for k in ['rawGateCount','rawGateCountLowerFullBins','rawGateCountUpperTouchingBins']}
previouswide=run.table(O/'atlas_low_002/candidate_atlas.tsv')
for response in cfg['responses']:
 folders={bg:R/('full_'+response+'_'+bg) for bg in cfg['backgrounds']};data={bg:run.read(d/'cutoff_fits.json')['slices'] for bg,d in folders.items()};warm={bg:run.read(R/('warm_'+response+'_'+bg+'/cutoff_fits.json'))['slices'] for bg in cfg['backgrounds'] if bg!='none'}
 settings[response]=run.read(folders['linear']/'effective_settings.json');settings[response]['fitWindowMode']='target_local';audits[response]={r['index']:r for r in run.read(folders['linear']/'tag_shape_audit.json')};group=collections.defaultdict(lambda:collections.defaultdict(list))
 for n,nobg in enumerate(data['none']):
  i=nobg['index'];variants={'none':(nobg,folders['none'],'anchor_start')}
  for bg in ['linear','exponential']:
   original,seeded=data[bg][n],warm[bg][n];options=[v for v in [(original,folders[bg],'anchor_start'),(seeded,R/('warm_'+response+'_'+bg),'none_start')] if valid(v[0])];variants[bg]=min(options,key=lambda x:x[0]['dataDeviance']) if options else (original,folders[bg],'failed')
  v=[x for x in previouswide if x['species']=='he4' and int(x['index'])==i];bounds={key:(min if key.endswith('Lo') else max)(float(x[key]) for x in v) for key in ['linearLo','linearHi','logLo','logHi']};bounds['linearLo']=min(bounds['linearLo'],nobg['fitLo']);bounds['logLo']=min(bounds['logLo'],nobg['fitLo'])
  for bg,(row,folder,start) in variants.items():
   model=response+'_'+bg;assess=old.aggregate_one([row],'he4',audits[response],settings[response],raw);c=row['components']['he4'];c3=row['components']['he3'];sel=c['selection'];record=dict(response=response,model=bg,index=i,r=row['r'],rLo=row['lo'],rHi=row['hi'],fitLo=row['fitLo'],fitHi=row['fitHi'],he3Mean=c3['mean'],he3Sigma=c3['sigma'],he4Mean=c['mean'],he4Sigma=c['sigma'],purity=sel['purity'],chi2=row['pearson'],ndf=row['conditionalNdf'],chi2Ndf=row['pearson']/row['conditionalNdf'] if valid(row) else None,deviance=row['dataDeviance'],numericalValid=valid(row),identityPass=row['targetIdentityScreenPass'],supported=assess['supported'],gateResidual=sel['localResidualFrac'],start=start,folder=str(folder),reasons=assess['reasons']);metrics.append(record)
   atlas.append(dict(band=response,species='he4',model=bg+'_'+start,index=i,rLo=row['lo'],rHi=row['hi'],fitLo=row['fitLo'],fitHi=row['fitHi'],mean=c['mean'],sigma=c['sigma'],root=str(folder/'cutoff_fits.root'),purity=sel['purity'],numericalValid=valid(row),supported=assess['supported'],**bounds,chi2=row['pearson'],ndf=row['conditionalNdf'],deviance=row['dataDeviance'],identityPass=row['targetIdentityScreenPass'],he3Mean=c3['mean'],he3Sigma=c3['sigma']))
   selected[model,i]=(row,folder);g=int(math.floor((row['lo']-cfg['aggregationOrigin'])/cfg['aggregateWidth']+1e-7));group[g][bg].append(row)
  models={bg:old.aggregate_one([r],'he4',audits[response],settings[response],raw) for bg,(r,p,start) in variants.items()}
  for bg in ['linear','exponential']:
   full=variants[bg][0];pair=dict(response=response,background=bg,index=i,r=row['r'],nestedViolation=valid(nobg) and valid(full) and full['dataDeviance']>nobg['dataDeviance']+cfg['nestedDevianceTolerance'],noneValid=valid(nobg),fullValid=valid(full),start=variants[bg][2]);details.append(pair)
 for g,models in sorted(group.items()):
  for bg,rows in models.items():aggregates.append(dict(response=response,model=bg,rLo=max(cfg['minimumRExclusive'],cfg['aggregationOrigin']+g*cfg['aggregateWidth']),rHi=min(cfg['tofThreshold'],cfg['aggregationOrigin']+(g+1)*cfg['aggregateWidth']),**old.aggregate_one(rows,'he4',audits[response],settings[response],raw)))
# Background-free fits are diagnostics. Final low-R purity has four nominal response x continuum models.
low=[]
for lo,hi in sorted({(r['rLo'],r['rHi']) for r in aggregates}):
 values={r['response']+'_'+r['model']:r for r in aggregates if r['rLo']==lo and r['rHi']==hi and r['model']!='none'};nominal=values['gaussian_linear'];supported=[v['purity'] for v in values.values() if v['supported'] and finite(v['purity'])];quality=len(values)==4 and all(v['supported'] for v in values.values());low.append(dict(mode='he4_r061',band='low',species='he4',rLo=lo,rHi=hi,purity=nominal['purity'],statConditionalBound=nominal['conditionalErrorWeightedBound'],purityLow=min(supported) if supported else None,purityHigh=max(supported) if supported else None,quality=quality,modelsSupported=sum(v['supported'] for v in values.values()),modelCount=4,status='supported_above_50' if quality and min(supported)>=cfg['purityThreshold'] else 'supported_below_50' if quality else 'unassessed',reasons=sorted({x for v in values.values() for x in v['reasons']}),models=values))
for i in sorted({i for model,i in selected}):
 models={model:old.aggregate_one([row],'he4',audits[model.split('_')[0]],settings[model.split('_')[0]],raw) for (model,j),(row,folder) in selected.items() if j==i and not model.endswith('_none')};rows=[row for (model,j),(row,folder) in selected.items() if j==i];supported=[v['purity'] for v in models.values() if v['supported'] and finite(v['purity'])];nominal=models['gaussian_linear'];native.append(dict(rLo=rows[0]['lo'],rHi=rows[0]['hi'],r=rows[0]['r'],purity=nominal['purity'],purityLow=min(supported) if supported else None,purityHigh=max(supported) if supported else None,quality=len(models)==4 and all(v['supported'] for v in models.values()),models=models))
prior=run.read(O/'candidate_summary.json');others=[r for r in prior['bins'] if not(r['species']=='he4' and r['band']=='low')];combined=others+low
run.save(R/'summary.json',dict(schemaVersion=1,productionEnabled=False,configuration=cfg,bins=combined,lowHe4=low,nativeHe4=native,metrics=metrics,backgroundComparisons=details,lowHe4Envelope='Four nominal Gaussian/lognormal x linear/exponential fits; no mean/width sensitivity variations. Other species retain their previous uncertainty envelope.',highHe4='Exact previous TOF result, with the deuteron m2 gate; no new TOF fits needed because PID gate values and high-R shapes are unchanged.'))
shared.tsv(R/'metrics.tsv',metrics,'response model index r rLo rHi fitLo fitHi he3Mean he3Sigma he4Mean he4Sigma purity chi2 ndf chi2Ndf deviance numericalValid identityPass supported gateResidual start folder'.split())
shared.tsv(R/'atlas.tsv',atlas,'band species model index rLo rHi fitLo fitHi mean sigma root purity numericalValid supported linearLo linearHi logLo logHi chi2 ndf deviance identityPass he3Mean he3Sigma'.split())
shared.tsv(R/'aggregate.tsv',aggregates,'response model rLo rHi purity signalSelected totalSelected rawSelected localResidualFrac supported'.split())
shared.tsv(R/'candidate_summary.tsv',combined,'mode band species rLo rHi purity statConditionalBound purityLow purityHigh quality modelsSupported modelCount status'.split())
shared.tsv(R/'he4_native_purity.tsv',native,'rLo rHi r purity purityLow purityHigh quality'.split())
for n in ['historical_summary.tsv','original_summary.tsv']:shutil.copy2(O/n,R/n)
shutil.copy2(O/'purity_axis_002/axes.tsv',R/'axes.tsv')
# Export the explicit study selection: same polynomial gate, new strict He4 minimum, corrected response calibration.
run.save(R/'pid_selection.yaml',dict(productionEnabled=False,coefficientSource=str(O/'coefficients.json'),coefficients=run.read(O/'coefficients.json')['candidate'],strictMinimumRigidityGeVcm={'he4':cfg['minimumRExclusive']},tofMass2=run.read(HERE.parent/'settings.yaml')['tofMass2'],responseCalibration=cfg))
print('SUMMARY',len(metrics),'native model rows;',sum(r['quality'] for r in low),'/',len(low),'low aggregate bins supported; nested violations',sum(r['nestedViolation'] for r in details))
for r in low[:3]:print('LOW_BIN',r['rLo'],r['rHi'],r['purity'],r['purityLow'],r['purityHigh'],r['quality'],r['reasons'])
