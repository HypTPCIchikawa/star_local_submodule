import datetime,fcntl,json,os,time
from pathlib import Path
S=Path(__file__).resolve().parent
lock=(S/'documentation.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
if (S/'documentation_completed.json').exists():raise SystemExit(0)
while True:
    try:
        result=json.loads((S/'result.json').read_text())
        if result['phase']=='failed' or result.get('archive',{}).get('phase') in ('archived','failed'):break
    except (OSError,ValueError):pass
    time.sleep(30)
for i in range(6):
    if (S/'sent.json').exists():break
    time.sleep(5)
repo=S.parents[4];tz=datetime.timezone(datetime.timedelta(hours=9));dt=datetime.datetime.now(tz)
good=result['phase']=='merged' and result.get('archive',{}).get('phase')=='archived'
label='完了（通常ROOT＋3 Sparse、個別ROOT保存済み）' if good else '未完了・要確認'
mail='ローカルsendmail受付済み' if (S/'sent.json').exists() else '送信受付を確認できていない'
text='\n\n## 自動完了記録 — '+dt.isoformat()+'\n\n状態：'+label+'。メール：'+mail+'。\n'
if result.get('merge'):text+='通常ROOT：`'+result['merge']['path']+'`（'+str(result['merge']['size'])+' bytes）。\n'
for part in result.get('sparseParts',[]):text+='Sparse：`'+part['path']+'`。\n'
if result.get('error'):text+='エラー：'+result['error']+'。\n'
text+='検証・公開・保存・送信記録：`'+str(S)+'`。\n'
public=repo/'rootfile'/('auau3p85_anaFemtoLambda_4He_KFParticle')/'README_outputs_04BF2DAE47164B6F040E3A460E539289.md'
with public.open('a') as out:out.write(text)
with (repo/'Todo/2026/10/26/he4_hadd_sparse_size_limit_20261026.md').open('a') as out:out.write(text)
todo_path=repo/'Todo/2026/10/26/STAR_Todo20261026.txt'
todo_text=todo_path.read_text()
if result['phase']=='merged':
    todo_text=todo_text.replace('[ ] 通常1,786ヒストの統合・最終検証・全4 ROOTの公開を完了する。','[x] 通常1,786ヒストの統合・最終検証・全4 ROOTの公開を完了する。')
if result.get('archive',{}).get('phase')=='archived':
    todo_text=todo_text.replace('[ ] 検証済み元26,955 ROOTをsubjobsへ保存する（削除しない）。','[x] 検証済み元26,955 ROOTをsubjobsへ保存する（削除しない）。')
if (S/'sent.json').exists():
    todo_text=todo_text.replace('[ ] 終了メールの送信受付と完了記録を確認する。','[x] 終了メールの送信受付と完了記録を確認する。')
todo_text=todo_text.replace('状態：β–η図・BOX重ね描きは完了。hadd/分割ROOTセットの最終処理は進行中。','状態：β–η図・BOX重ね描きは完了。hadd/分割ROOTセット：'+label+'。')
if good:todo_text=todo_text.replace('■ 1. 通常ヒスト統合＋巨大PID Sparseの3分割保存【進行中】','■ 1. 通常ヒスト統合＋巨大PID Sparseの3分割保存【完了】')
todo_path.write_text(todo_text+text)

log=repo/'analysisnote'/dt.strftime('%Y%m%d');log.mkdir(exist_ok=True);summary=log/('summary'+dt.strftime('%Y%m%d')+'.md')
if not summary.exists():summary.write_text('# 解析ノート '+dt.strftime('%Y-%m-%d')+'\n')
with summary.open('a') as out:out.write('\n\n---\n\n## Session update — '+dt.strftime('%H:%M JST / %Y-%m-%d')+'：分割出力の完了記録\n'+text)
(S/'documentation_completed.json').write_text(json.dumps({'writtenAt':dt.isoformat(),'phase':result['phase'],'archive':result.get('archive',{}).get('phase'),'mailAccepted':(S/'sent.json').exists()},indent=2)+'\n')
