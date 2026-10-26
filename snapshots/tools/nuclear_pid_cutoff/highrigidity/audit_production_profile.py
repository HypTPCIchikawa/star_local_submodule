#!/usr/bin/env python3
"""Independent production-profile/rollback audit; source/config are read-only."""
import argparse
import hashlib
import json
import shutil
from pathlib import Path
import yaml
import run_highrigidity as high
common=high.common


def accepted(minimum,veto,maximum):
    result=[];edge=minimum
    for interval in veto:
        lo,hi=interval['min'],interval['max']
        if lo<edge or hi<=lo or hi>maximum: raise ValueError('Unordered/overlapping/out-of-domain veto.')
        if lo>edge: result.append([edge,lo])
        edge=hi
    if edge<maximum: result.append([edge,maximum])
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--domains',required=True,type=Path)
    parser.add_argument('--profile',required=True,type=Path)
    parser.add_argument('--before-dir',required=True,type=Path)
    parser.add_argument('--reference-run',required=True,type=Path)
    parser.add_argument('--mains',required=True,nargs='+',type=Path)
    parser.add_argument('--rollback-suffix',required=True)
    parser.add_argument('--output-dir',required=True,type=Path)
    args=parser.parse_args();domains_dir=common.resolve(args.domains);profile=common.resolve(args.profile)
    before=common.resolve(args.before_dir);reference=common.resolve(args.reference_run)
    done=common.read(domains_dir/'completed.json');domains=common.read(domains_dir/'production_domains.json')
    for name,sha in done['products'].items():
        if common.digest(domains_dir/name)!=sha: raise ValueError('Changed domain artifact: '+name)
    for path,sha in domains['protectedInputs'].items():
        if common.digest(Path(path))!=sha: raise ValueError('Changed immutable range input: '+path)
    cuts=yaml.safe_load(profile.read_text());old=common.read(reference/'current_selection.json')
    old_context=common.read(reference/'context.json')
    if set(cuts)-set(old['cuts']): raise ValueError('Unmodeled new production keys.')
    changed={k for k in cuts if cuts[k]!=old['cuts'].get(k)}
    if changed-{'maxPOverQ','speciesPOverQVeto','highRigidityTof'}: raise ValueError('Unexpected non-domain PID change: '+repr(changed))
    if cuts['maxPOverQ']!=domains['maximum'] or cuts['m2_selection'] or cuts['requireBestSpecies'] is not True:
        raise ValueError('Original TPC-winner / post-selection policy changed.')
    if set(cuts['speciesMinPOverQ'])!=set(high.SPECIES) or set(cuts['speciesPOverQVeto'])!=set(high.SPECIES):
        raise ValueError('Missing/unknown species domain.')
    checks={}
    for species in high.SPECIES:
        expect=domains['domains'][species]
        got=accepted(cuts['speciesMinPOverQ'][species],cuts['speciesPOverQVeto'][species],cuts['maxPOverQ'])
        if cuts['speciesMinPOverQ'][species]!=expect['minimum'] or cuts['speciesPOverQVeto'][species]!=expect['veto'] or got!=expect['accepted']:
            raise ValueError('Production YAML differs from accepted study domain: '+species)
        checks[species]=got
    tof=cuts['highRigidityTof'];old_tof=old['cuts']['highRigidityTof']
    if tof['enabled'] is not True or tof['thresholdPOverQ']!=domains['split'] or tof['extendHe4BeyondMaxPOverQ'] is not False:
        raise ValueError('Unexpected TOF threshold/global-domain policy.')
    expected_tof=dict(old_tof['species'],he4=old_tof['species']['d'])
    if tof['species']!=expected_tof: raise ValueError('TOF windows differ from fitted exact samples.')
    sources={}
    for path,sha in old['protectedFiles'].items():
        if path in (old['runtimeMainconf'],old['nuclearConfig']): continue
        sources[path]=dict(unchanged=common.digest(Path(path))==sha,expectedSha256=sha)
    for path,sha in old_context['protectedSources'].items():
        if path in (str(common.PROJECT/'StMaker/common/FemtoLambdaLegacy.cxx'),str(common.PROJECT/'StMaker/common/FemtoLambdaLegacy.h')):
            sources[path]=dict(unchanged=common.digest(Path(path))==sha,expectedSha256=sha)
    if not all(v['unchanged'] for v in sources.values()): raise ValueError('TPC provider/calibration or winner implementation changed.')
    if len(args.mains)!=5 or len({p.name for p in args.mains})!=5: raise ValueError('Exactly five distinct live mains required.')
    main_checks={};protected={str(profile):common.digest(profile)}
    expected_ref=str(profile.relative_to(common.PROJECT/'config'))
    for supplied in args.mains:
        main=common.resolve(supplied);previous=before/main.name
        rollback=main.with_name(main.stem+args.rollback_suffix+main.suffix)
        a,b,c=(yaml.safe_load(p.read_text()) for p in (main,previous,rollback))
        if a.get('nuclearid')!=expected_ref: raise ValueError('Live main does not use new shared profile: '+main.name)
        if {k:v for k,v in a.items() if k!='nuclearid'}!={k:v for k,v in b.items() if k!='nuclearid'}:
            raise ValueError('Unexpected non-nuclearid main change: '+main.name)
        if b!=c: raise ValueError('Rollback is not semantically identical to before: '+rollback.name)
        for key,value in a.items():
            if not isinstance(value,str): raise ValueError('Unsupported main reference: '+key)
            linked=(common.PROJECT/'config'/value).resolve(strict=True)
            protected[str(linked)]=common.digest(linked)
        analysis=yaml.safe_load((common.PROJECT/'config'/a['analysis']).read_text())
        tag=analysis['starTag']
        if tuple(tag[k] for k in ('collision','year','productionTag','libraryTag'))!=('auau3p85','year2021','P24iy','SL24y'):
            raise ValueError('A main points outside 2021 P24iy/SL24y.')
        for p in (main,previous,rollback): protected[str(p)]=common.digest(p)
        main_checks[main.name]=dict(before=str(previous),rollback=str(rollback),nuclearidOnlyChange=True,
            rollbackSemanticallyIdentical=True,allReferencesExist=True,dataset2021P24iy=True,
            otherReferences={k:v for k,v in a.items() if k!='nuclearid'},previousNuclearid=b['nuclearid'],analysisName=analysis['analysis']['anaName'])
    before_cuts={}
    for old_path in before.glob('nuclearid*.yaml'):
        original=common.PROJECT/'config/cuts/nuclearid'/old_path.name
        same=yaml.safe_load(old_path.read_text())==yaml.safe_load(original.read_text())
        if not same: raise ValueError('Old rollback nuclear settings changed: '+str(original))
        before_cuts[old_path.name]=dict(unchanged=True,sha256=common.digest(original))
    output=args.output_dir.resolve();output.mkdir(parents=True,exist_ok=True)
    if (output/'completed.json').exists() or (output/'audit.json').exists(): raise ValueError('Existing audit will not be overwritten.')
    report=dict(schemaVersion=1,passed=True,acceptedDomains=checks,physicalProfileChangedKeys=sorted(changed),
        beforeNuclearSettingsUnchanged=before_cuts,currentMainChecks=main_checks,providerAndWinnerSourcesUnchanged=sources,
        highTofWindows=tof['species'],postWinnerDomainOrderPreserved=True,
        dAndTCandidatesMayChangeBecauseOfCompetingHeliumCalibration=True,
        he3OwnResponseMayChangeFromLegacyToSharedCalibratedProvider=True,
        expectedSameYieldAsOldAnalyses=False,note='Smoke validation checks configuration and execution, not equal physics yields or efficiency.')
    high.save(output/'audit.json',report)
    high.save(output/'source_inputs.json',protected)
    source=output/'source';source.mkdir();copy=source/Path(__file__).name
    shutil.copy2(Path(__file__),copy)
    for path,sha in protected.items():
        if common.digest(Path(path))!=sha: raise ValueError('Input changed during audit.')
    high.save(output/'completed.json',dict(executionComplete=True,independentAuditPassed=True,
        products={str(p.relative_to(output)):common.digest(p) for p in output.rglob('*') if p.is_file()}))
    print('PRODUCTION_PROFILE_INDEPENDENT_AUDIT_PASS mains=5 rollback=5 domains=4',flush=True)
