#!/usr/bin/env python3
import json,pathlib,re,shutil,subprocess,sys
HERE=pathlib.Path(__file__).resolve().parent;sys.path.insert(0,str(HERE));import prepare as base
R=base.R;C=base.C;summary=base.read(R/'summary.json');coef=base.read(R/'coefficients.json');examples=base.read(R/'examples.json')
def text(s):
 replacements={'\\':r'\textbackslash{}','&':r'\&','%':r'\%','$':r'\$','#':r'\#','_':r'\_','{':r'\{','}':r'\}','~':r'\textasciitilde{}','^':r'\textasciicircum{}'}
 s=''.join(replacements.get(c,c) for c in s)
 s=s.replace('/',r'/\allowbreak{}').replace(r'\_',r'\_\allowbreak{}')
 return re.sub(r'([\u2e80-\u9fff\uf900-\ufaff\uff00-\uffef]+)',lambda m:r'\hspace{0pt}{\jpfont '+m.group(1)+r'}\hspace{0pt}',s)
def p(s):return text(s)+'\n\n'
def head(s):return r'\section*{'+text(s)+'}\n'
def fig(name,height='.60'):
 name=name.replace('.png','_trim.png')
 return r'\begin{center}\includegraphics[width=\linewidth,height='+height+r'\textheight,keepaspectratio]{'+name+r'}\end{center}'+'\n'
def table(headers,rows):
 align='l'+'r'*(len(headers)-1);lines=[r'\begin{center}\small\begin{tabular}{'+align+r'}\toprule',' & '.join(text(v) for v in headers)+r'\\\midrule']
 lines+=['{}'+' & '.join(text(str(v)) for v in row)+r'\\' for row in rows];return '\n'.join(lines)+r'\bottomrule\end{tabular}\end{center}'+'\n'
first=next(r for r in summary['rows'] if r['species']=='he4' and r['band']=='low');pages=[]
pages.append(head('原子核PIDとPurity推定')+p('Gaussian応答＋linear背景による最終表示 / 2026年10月21日')+p('対象：2021 P24iy、auau3p85解析試料。d、t、³He、⁴Heについて、R=|p/q|ごとのPurityを推定する。R<2 GeV/cはTOFを使わないALL、R≥2 GeV/cは既存TOF m²選択を用いる。')+p('目視判断とユーザー指定に従い、全核種をGaussian＋linearの名目fitだけに統一した。背景なし・lognormal・exponentialの点、他モデル由来の範囲はこのPurity図に含めない。')+fig('purity_vs_p_over_q_01.png','.57')+p('図1　更新した最終Purity。青：ALL、赤：TOF。横誤差はR区間、縦誤差は収量の条件付き誤差。灰色の×はGaussian＋linear自身の支持条件を通らない区間。')+p('⁴HeはR>0.61を適用。0.61–0.65のPurityは'+f"{100*first['purity']:.2f} ± {100*first['statConditionalBound']:.2f}"+'%（条件付き誤差）。本報告は解析候補の更新をまとめたもので、本番設定への反映は行っていない。'))
pages.append(head('1. 選別とPurityの定義')+p('Rは電荷で割った運動量であり、電荷2のHeでは運動量そのものと異なる。TPCのdE/dx PIDは中心μPID(R)と幅σPID(R)で定義し、厳密な|nσ|<2とbest-species判定を用いる。重複時は|nσ|が最小の核種を採用し、同値の優先順はd、t、³He、⁴He。')+r'\[n\sigma_i(R,x)=\frac{x-\mu_{\mathrm{PID},i}(R)}{\sigma_{\mathrm{PID},i}(R)},\qquad |n\sigma_i|<2.\]'+'\n'+p('観測スペクトルはd、t、³He、⁴HeのGaussian成分とlinearなsmooth背景の和でfitする。linear背景はfit領域内で非負になる端点基底を用いる。右側のlog表示は縦軸の表示方法であり、lognormal応答を採用したものではない。')+r'\[f(x;R)=\sum_{i=d,t,{}^3\mathrm{He},{}^4\mathrm{He}}N_i(R)G_i(x;\mu_i,s_i)+B_{\rm linear}(x;R).\]'+'\n'+p('核種iのPID選択領域へ各成分を積分し、Sを目的核種、Bを他核種の混入とlinear背景の和とする。複数native区間をまとめる際はPurityの平均ではなく、分子・分母を先に合計する。')+r'\[P_i=\frac{S_i}{S_i+B_i},\qquad P_i(\Delta R)=\frac{\sum_k S_{i,k}}{\sum_k (S_{i,k}+B_{i,k})}.\]'+'\n'+table(['R領域','試料・条件'],[['R<2','ALL：TOFによる試料選択なし'],['R≥2','既存の核種別TOF m²窓'],['⁴He全域','strict R>0.61、TOF側はdと同じm²窓']])+table(['核種','TOF m²下限','TOF m²上限'],[[s,*[f'{v:.7g}' for v in pair]] for s,pair in base.read(HERE.parent/'settings.yaml')['tofMass2'].items()])+p('通常の表示は0.05 GeV/c幅。低Rのnative幅は0.01。⁴Heの最初の集約区間のみ下限でclipし、0.61–0.65の4 native区間を使う。TOF側⁴Heは既存の0.05幅fitを維持する。R方向はGL3で積分し、選択窓の生計数はnative dE/dx binとの重なり比で評価する。イベント単位のPID再選択ではない。'))
meanrows=[];sigrows=[]
for i,s in enumerate(['d','t','³He','⁴He']):
 meanrows.append([s]+[f'{x:.9g}' for x in coef['candidate']['mean'][i]]);sigrows.append([s]+[f'{x:.9g}' for x in coef['candidate']['sigma'][i][:3]])
pages.append(head('2. PID関数形と係数の更新')+p('Git 8f7d539の元の関数形へ戻り、x=log10 Rに対する中心3次・幅2次を維持した。d/tの係数は不変更。Heの係数のみ既存ALL・TOFの自由fit参照に基づき微調整した。')+r'\begin{align*}x&=\log_{10}\!\left(R/(1\,{\rm GeV}/c)\right),\\ \mu_{\rm PID}(R)&=a_0+a_1x+a_2x^2+a_3x^3,\\ \sigma_{\rm PID}(R)&=b_0+b_1x+b_2x^2.\end{align*}'+'\n'+table(['核種','a0','a1','a2','a3'],meanrows)+table(['核種','b0','b1','b2'],sigrows)+p('係数の調整では中心変化を元のσの0.5倍以内、幅変化を元のσの20%以内とし、R=0.4–5の監査域で中心と幅の正値を検証した。表は丸めた表示であり、再現には同梱coefficients.jsonの全桁を用いる。')+p('重要な区別：このμPID、σPIDは選択窓の関数であり、各1次元スペクトルfitのGaussian中心μfit、幅sfitと同一ではない。後述の低Rピーク修正はfitの参照成分に反映しており、TPC PIDの多項式を局所補間へ置き換えていない。')+p('今回追加した選択条件は⁴HeのR>0.61。多項式の中心・幅の係数は直前の候補から維持した。³Heの選択関数を、⁴He fit内の³He参照中心と混同しない。')+p('本番PID設定は未変更。解析候補の選択条件、係数元と低R応答較正をpid_selection.yamlに分けて保存している。'))
pages.append(head('3. PID中心の変更前後')+fig('pid_function_updates_01.png','.72')+p('図2　元のGit関数（灰色破線）と解析候補（色線）の中心μPID。d/tは重なる。Heは従来の3次式の係数だけを微調整した。⁴Heの候補曲線は適用範囲R>0.61に限定して表示した。')+p('この図はPID選択窓の中心であり、低R Gaussian fitの中心を示す図ではない。指定されたR=0.61/0.64のピーク位置は後述の応答較正へ反映する。'))
pages.append(head('4. PID幅の変更前後')+fig('pid_function_updates_02.png','.72')+p('図3　元のGit関数と候補の幅σPID。幅はlog10 Rの2次式であり、新しいガタガタした補間は導入していない。幅の正値と変更量を監査した。')+p('選択窓はμPID±2σPIDに加え、核種間のbest-species判定で区切られる。したがって、図の帯を単独で切り取っただけの選択ではない。'))
pages.append(head('5. dE/dx分布と選択曲線')+fig('dedx_pid_map_01.png','.61')+p('図4　ALLのdE/dx対R。灰色はGit元の曲線、色線は更新候補の中心と±2σ。低R側のPurity計算はこのALL試料で行い、TOFで試料を狭めない。⁴He候補はR>0.61を適用する。')+p('曲線の表示とfit範囲は別である。1次元fitの横軸は両ピークを観察できるよう拡張し、実際のfit下限は³Heピークの低側半値を含むよう設定した。fit範囲外へ予測曲線を外挿していない。')+p('PID多項式の調整と、ALL中のGaussian応答中心・幅の再fitを分けて扱うことで、選択窓の変更と応答仮定の影響を追跡できる。'))
anch=[r for r in examples if r['group']=='anchors'];arows=[[f"[{r['rLo']:.2f},{r['rHi']:.2f})",f"{r['he3Mean']:.2f}",f"{r['mean']:.2f}",f"{100*r['purity']:.2f}%"] for r in anch]
pages.append(head('6. ⁴He低Rピークの修正')+p('指定参照はR=0.61で³He=51、⁴He=67、R=0.64で³He=49、⁴He=65 keV/cm。³He中心はこの較正参照へ固定し、0.61–0.70だけ幅をALLから自由fitした。⁴He中心と幅はALLから自由fitし、指定値へ強制固定していない。')+fig('fit_anchors_01.png','.52')+p('図5　上：[0.61,0.62)、下：[0.64,0.65)。左右はlinear/log表示で、どちらもGaussian＋linear fit。native中心はそれぞれ0.615、0.645であり、指定Rの無幅分布とは異なる。')+table(['R区間','³He中心・固定','⁴He中心・自由','Purity'],arows)+p('³Heの低側半値を全fit・全GL3点で含むことを、fit後の幅についても検証した。0.70以上の1,300比較は直前の同じ半値範囲fitと数値差ゼロ。'))
pages.append(head('7. 低Rの代表fit：dとt')+fig('fit_low_01.png','.68')+p('図6　d、tの代表例。左linear、右log。黒は生データ、赤は全成分、灰破線はlinear背景、色線は核種成分。すべて保存済みGaussian＋linearのfitから再描画し、新しいfitは実行していない。')+p('数値や曲線はfit領域内だけに表示する。図のχ²/NDF、Purity、中心・幅から、形状と選択領域での混入を併せて確認できる。'))
pages.append(head('8. 低Rの代表fit：³Heと⁴He')+fig('fit_low_02.png','.68')+p('図7　³He、⁴Heの代表例。核種のGaussian成分とlinear背景を分けて示す。⁴He fit中の³He成分は混入評価に寄与するため、⁴Heピークだけが見える狭い表示は避けた。')+p('R<2ではTOFを使わずALLで評価する。linear背景を選んだことは、他核種の成分を背景計算から除く意味ではない。Purityの分母には選択窓へ入る全成分を含める。'))
pages.append(head('9. TOF領域の代表fit：dとt')+fig('fit_tof_01.png','.68')+p('図8　R≥2のTOF選択試料。既存m²窓による選択とGaussian＋linearのfit結果をそのまま使った。低RのALLと高RのTOFは試料選択が異なるため、R=2の境界を明示してPurityを接続する。')+p('高Rの支持判定には、従来の数値収束、応答誤差、収量、核種成分を除いた際のdeviance差等を用いる。未評価点を補間して延長しない。'))
pages.append(head('10. TOF領域の代表fit：³Heと⁴He')+fig('fit_tof_02.png','.68')+p('図9　³He、⁴HeのTOF領域。⁴Heはdとm²が重なるため、dと同一のTOF選択窓を使用する。その後にdE/dxの成分fitとPID積分で⁴Heを評価する。')+p('TOF m²だけでdと⁴Heを分離できたと扱わず、両核種のGaussian成分を含めてPurityを計算する。高Rの元のfit・入力ROOTは変更していない。'))
highstats=[r for r in summary['statistics'] if r['band']=='high'];ht=[[r['species'],f"{r['supported']}/{r['total']}",f"{r['rMin']:.2f}–{r['rMax']:.2f}",f"{100*r['minPurity']:.2f}%"] for r in highstats]
pages.append(head('11. 高RのPurity')+fig('purity_vs_p_over_q_02.png','.56')+p('図10　Gaussian＋linearのみの高R表示。誤差棒は条件付き収量誤差で、モデル差の範囲ではない。')+table(['核種','支持点数/全区間','支持点のR最小–最大','最小Purity'],ht)+p('表示範囲は各核種の支持点に合わせた。dの全体図の横軸は従来の0.4–5を保持。支持点の最大Rはこの解析での推定可能範囲であり、物理的・統計的な絶対上限の確定ではない。'))
below=[[r['species'],f"{r['rLo']:.2f}–{r['rHi']:.2f}",f"{100*r['purity']:.2f}%"] for r in summary['rows'] if r['quality'] and r['purity']<C['purityThreshold']]
pages.append(head('12. 採否、誤差と再現性')+p('最終図の採否はGaussian＋linear自身の支持条件だけを使う。以前の複数モデル全体の採否や包絡線を残さず、他モデルの失敗を理由に採用モデルの点を落とさない。一方、Purityが50%未満の点は隠さず表示する。')+table(['核種','集約R区間','Gaussian＋linear Purity'],below)+p('⁴HeのR>0.61は適用範囲を表し、全Rで50%以上を保証するものではない。特に1.95–2.00では45.04%が残る。³Heの低Rにも50%未満の区間がある。新しいvetoは自動追加していない。')+p('縦誤差はfit形状・背景モデルを条件とした収量誤差。集約では各native区間の誤差を分母の重みで線形に足した保守的な上限を用いる。Gaussian/linearの選択、応答形状、³He参照位置の概数誤差を含む完全な系統誤差ではない。単一モデル表示に変更して誤差棒が小さくなっても、全不確かさが減ったという意味ではない。')+p('絶対χ²/NDFが大きい区間もある。支持フラグの通過はχ²/NDF≈1を意味しない。linear採用はユーザーによる目視判断と指定を反映したもので、全域での形式的なモデル優越を主張していない。')+p('今回の作業は既存fitの抽出と再描画で、新しいfit・イベント解析・本番PID反映は行っていない。364区間のモデルキー、ΣS/Σtotal、条件付き誤差と支持フラグを元JSONへ照合する。旧ROOT・旧図を保持し、入力SHAと生成物をmanifestに記録する。')+p('再現ファイル：同梱settings.yaml、purity.tsv、summary.json、coefficients.json、examples.tsv。ツールはtools/nuclear_pid_smooth_20261021/linear_report/。元はcandidate_summary.jsonとhe4_r061_001/summary.json。')+p('mainconf：main_auau3p85_anaFemtoLambda_4He_KFParticle.yaml。描画はscript/singularity_root_batch.sh経由のSTAR SL24y / ROOT5。ローカルstudyのため新規jobid/configlogなし。PDF組版はXeLaTeX。図のvector PDFと200 dpi画像を両方保存。job投入・hadd・commit/pushなし。'))
preamble=r'''\documentclass[10pt,a4paper]{article}
\usepackage[margin=18mm,top=20mm,bottom=18mm]{geometry}
\usepackage{fontspec,graphicx,amsmath,booktabs,xcolor,fancyhdr}
\setmainfont{DejaVu Sans}
\newfontfamily\jpfont[Path=/usr/share/fonts/google-droid-sans-fonts/,BoldFont=DroidSansFallbackFull.ttf,BoldFeatures={FakeBold=1.2}]{DroidSansFallbackFull.ttf}
\XeTeXlinebreaklocale "ja"
\XeTeXlinebreakskip=0pt plus 1pt
\definecolor{navy}{HTML}{17385D}
\pagestyle{fancy}\fancyhf{}\fancyhead[L]{\small STAR / Nuclear PID and Purity}\fancyhead[R]{\small Gaussian + linear}\fancyfoot[C]{\thepage}\setlength{\headheight}{14pt}
\emergencystretch=2em\setlength{\parindent}{0pt}\setlength{\parskip}{5pt}\linespread{1.10}
\begin{document}
'''
body='\n\\clearpage\n'.join(pages);tex=preamble+body+'\n\\end{document}\n';(R/'pid_purity_report_20261021.tex').write_text(tex)
for attempt in range(2):
 cmd=['xelatex','-interaction=nonstopmode','-halt-on-error','-file-line-error','pid_purity_report_20261021.tex']
 with (R/('latex_build_'+str(attempt+1)+'.log')).open('w') as f:status=subprocess.run(cmd,cwd=R,stdout=f,stderr=subprocess.STDOUT).returncode
 assert not status,'XeLaTeX failed; inspect build log'
log=(R/'pid_purity_report_20261021.log').read_text();assert 'Missing character:' not in log,'Missing glyphs';base.save(R/'report_build.json',dict(expectedPages=len(pages),language='Japanese',figurePages=10,figureDpi=C['figureRasterDpi'],command=cmd,missingGlyphs=False,overfullBoxes=log.count('Overfull \\hbox')))
print('REPORT_BUILT',len(pages),'designed pages; overfull boxes',log.count('Overfull \\hbox'),flush=True)
