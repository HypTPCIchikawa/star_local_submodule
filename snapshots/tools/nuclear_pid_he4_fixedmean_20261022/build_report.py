#!/usr/bin/env python3
import pathlib,sys,subprocess,re
HERE=pathlib.Path(__file__).resolve().parent;sys.path.insert(0,str(HERE));import run
R,C=run.R,run.C;summary=run.read(R/'summary.json');coef=run.read(R/'coefficients.json');native=run.read(R/'native.json');validation=run.read(R/'physics_validation.json');STEM='he4_pid_purity_report_20261022'
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
pages=[];first=summary['rows'][0];stats=summary['statistics'];valid=[r for r in summary['rows'] if r['quality']];low=[r for r in valid if r['band']=='low'];high=[r for r in valid if r['band']=='high']
pages.append(head('⁴He Purity再評価：低R中心固定')+p('2026年10月22日 / 2021 P24iy、auau3p85解析試料')+p('改訂版：R<0.8のlinear背景を除き、前版Gaussian＋linearの⁴He中心に固定して再評価。PID係数とR≥0.8の結果は保持した。')+p('既存の2021 P24iy部分統合試料を使い、対象を⁴Heだけに絞る。d、t、³HeのPID関数形・係数・元の適用域はGit 8f7d539に固定し、⁴Heだけ従来多項式の係数を調整した。低RとTOF領域を指定の異なるモデルで評価する。')+fig('purity_vs_p_over_q_01.png','.56')+p('図1　⁴Heの最終Purity。青：ALLのGaussianのみ・中心固定、緑：ALLのGaussianのみ・中心自由、赤：TOFの⁴He Gaussian＋linearのみ。灰色×は未支持区間の位置を下端へ示し、Purity=0を意味しない。境界で制約・モデル・試料が変わるため、点を線で接続していない。') +p(f"R=0.61–0.65のモデル計算値は{100*first['purity']:.5f}%だが、幅の支持条件を通らないため採用点に含めない。全88集約区間中{len(valid)}区間が支持条件を通過。本番設定への反映は行っていない。"))
pages.append(head('1. 選択と3領域の評価モデル')+table(['R [GeV/c]','試料','fitの成分'],[['0.61–0.80','ALL、TOFなし','4核種Gaussianのみ・⁴He中心固定'],['0.80–2.00','ALL、TOFなし','d/t/³He/⁴He Gaussianのみ'],['2.00–5.00','既存TAG_D','⁴He Gaussian＋linearのみ']])+p('R=|p/q|。⁴Heはstrict R>0.61、0.80は⁴He中心自由側に含める。TOF側はR=2.00からのnative区間を使う。境界点R=2自体は連続分布上で測度0であり、区間に重複・欠落はない。⁴HeのTOF m²窓はdと同一のstrict 3.268773<m²<3.693147。')+p('TPC選択は全4核種の|nσ|<2とbest-species判定。競合時は|nσ|最小、同値の優先順d、t、³He、⁴He。d/t/³He自身の式が固定でも、⁴Heの係数変更により競合境界は変わるため、新しい全核種の選択窓を使って積分し直した。')+r'\[n\sigma_i=\frac{dE/dx-\mu_{{\rm PID},i}(R)}{\sigma_{{\rm PID},i}(R)},\qquad P_4=\frac{S_4}{S_4+B}.\]'+'\n'+p('低RのBは選択窓に入る他核種のGaussian寄与。R<2にsmooth背景はない。高RのBはlinearだけであり、他核種のGaussian成分はfitにも分母にも入れない。これはユーザー指定モデルの仮定であり、データから他核種が不存在と証明したものではない。')+r'\[\text{TOF:}\quad f=N_4G_4+B_{\rm linear},\qquad P_4=\frac{\int_{\rm PID}N_4G_4}{\int_{\rm PID}(N_4G_4+B_{\rm linear})}.\]'+'\n'+p('低R native幅0.01、TOF側は既存0.05 GeV/c。Purity表示は0.05幅へΣS/Σtotalで集約し、最初のみ0.61–0.65へclipする。R方向はGL3で積分、生計数はnative dE/dx binの重なり比から独立計数する。'))
mu=[[label]+[f'{v:.9g}' for v in coef['original']['mean'][i]] for i,label in enumerate(['d・固定','t・固定','³He・固定','⁴He・元'])];mu.append(['⁴He・候補']+[f'{v:.9g}' for v in coef['candidate']['mean'][3]])
sigma=[[label]+[f'{v:.9g}' for v in coef['original']['sigma'][i][:3]] for i,label in enumerate(['d・固定','t・固定','³He・固定','⁴He・元'])];sigma.append(['⁴He・候補']+[f'{v:.9g}' for v in coef['candidate']['sigma'][3][:3]])
pages.append(head('2. Git固定値と⁴He候補の係数')+r'\[x=\log_{10}(R/(1\,{\rm GeV}/c)),\quad\mu_{\rm PID}=a_0+a_1x+a_2x^2+a_3x^3,\quad\sigma_{\rm PID}=b_0+b_1x+b_2x^2.\]'+'\n'+table(['核種','a0','a1','a2','a3'],mu)+table(['核種','b0','b1','b2'],sigma)+p('Gitヘッダを直接保存し、d/t/³Heの係数・関数形・適用域の一致を検証した。以前のstudyで微調整した³He係数は使用していない。候補の全桁はcoefficients.jsonとcandidate_pid/NuclearIdDeDxVsMom.hに保存し、Gitとの差分も同梱する。')+p('前版の基準3領域fitの199点中169点を、数値成立・形状支持・局所残差とR=0.61–3.5で選んで係数を推定。元のσ単位で等しい重みの最小二乗を用い、Purity最大化による調整はしない。中心は元σの0.5倍、幅は元σの20%以内へ変化を制限した。')+p('今回、PID係数を再推定せず全桁保持した。前版の監査域R=0.61–5の1,001点で正の中心・幅と変更量を確認。中心の正規化RMSは0.268→0.210、幅は0.218→0.146。PID関数とスペクトル応答の中心・幅は別の量である。'))
pages.append(head('3. ⁠⁴Heの中心・幅の変更前後')+fig('pid_functions_01.png','.71')+p('図2　灰色破線はGit元PID、紫は⁴He候補。上段は全R、下段は低R。青の点は支持された各スペクトルfitのGaussian応答で、PID選択窓の中心・幅と同一視しない。d/t/³HeはGit元係数のまま。')+p('曲線は元の多項式で滑らかに定義し、R=0.8や2.0でPID関数を切り替えていない。今回0.80で変わるのは⁴He応答中心の固定/自由、2.00で変わるのはTOF要求とlinear背景・核種成分である。'))
pages.append(head('4. ALL分布とPID選択曲線')+fig('dedx_pid_map_01.png','.62')+p('図3　TOFを要求しないALLのdE/dx対R。d/t/³HeはGit元の中心と±2σ、⁴Heは灰色の元関数と紫の候補を示す。実際の選択は帯の単独切出しに加えbest-speciesで区切る。')+p('高Rで他核種のfit成分を除いたことは、TPC PIDの競合判定からd/t/³Heを取り除いたことを意味しない。元の競合判定を保った⁴He選択に対し、指定モデルでPurityを求めた。')+p('⁴Heの適用下限R>0.61はpid_selection.yamlで明示する。候補ヘッダの元の関数適用域自体は保持しており、ヘッダだけで下限cutまで実装したものではない。'))
pages.append(head('4a. TOF m²選択後の分布とPID選択曲線')+fig('dedx_tof_pid_map_01.png','.68')+p('図3a　既存TAG_D試料のR–dE/dx分布。TOFのstrict 3.268773<(m/q)²<3.693147を要求し、TPC nσ選択・best-species判定の前の分布を使った。上段はR=0.61–5、下段は高R=2–5 GeV/cの拡大。色は対数スケールのtrack数。')+p('中心と±2σのPID曲線を重ねた。d/t/³HeはGit元係数、⁴Heは灰色破線がGit元、紫が今回確定した係数。本reportの既存候補曲線と同じ全桁のcoefficients.json/curves.tsvを使用した。高RのPurity評価で用いたTOF試料の核バンドとPID位置を確認するための図であり、データの再fitはしていない。')+p('低RにもTOF選択済み分布を参考表示するが、R<2の最終Purityは従来どおりTOF不要のALLで評価する。既存2021 P24iy部分統合試料由来のcacheを使い、現在実行中の本番job出力は使用していない。'))
for sec,region in enumerate(C['regions'],5):
 name=region['name'];title={'low_fixed':'低R Gaussianのみ・中心固定','low_none':'中R Gaussianのみ','high_linear':'高R ⁠⁴He Gaussian＋linearのみ'}[name]
 reps=[r for r in run.read(R/'representatives.json') if r['region']==name];intervals='、'.join(f"[{r['rLo']:.2f},{r['rHi']:.2f})" for r in reps)
 intro='⁴He中心は前版Gaussian＋linearの各R区間の値へ固定し、幅と核種収量は自由fit。低R³Heの較正条件も保持し、linear背景を完全に除去した。' if name=='low_fixed' else 'smooth背景を完全に除く。Gaussianは他核種の混入も含む4成分であり、⁴He信号だけのfitではない。' if name=='low_none' else 'TOFはdと同じm²窓。fitは⁴Heとlinear背景だけで、d/t/³Heの寄与は0、自由パラメータにも含めない。'
 pages.append(head(f'{sec}. 代表fit：{title}')+p('R='+intervals+' GeV/c。各段で左linear、右log表示。'+intro)+fig('fit_representative_'+name+'_01.png','.68')+p(f'図{sec-1}　黒：生データ、赤：全モデル、灰破線：linear背景（ある領域のみ）、色線：核種成分。χ²/NDF、Purity、fit範囲、中心・幅、数値成立と支持フラグを表示。予測曲線はfit範囲内のみ描画した。')+p('低Rではfit後の³He低側半値を全GL3点で含む。右のlogは縦軸表示で、lognormal応答を用いたものではない。' if name!='high_linear' else '高Rは前版で⁴Heとlinearだけの60区間を再fitしており、今回はその結果を全桁保持した。非表示で隠しただけの核種成分はない。'))
st=[[{'low_fixed':'ALL中心固定・背景なし','low_none':'ALL背景なし','high_linear':'TOF ⁠⁴He＋linear'}[r['region']],f"{r['supported']}/{r['total']}",f"{100*r['minPurity']:.2f}%"] for r in stats]
pages.append(head('8. 低R Purityと評価可能域')+fig('purity_vs_p_over_q_02.png','.57')+p('図7　低R表示。0.80で⁴He応答中心の固定/自由を切り替える。0.61–2.00は全域でsmooth背景なし。')+table(['領域','支持/全区間','支持点の最小Purity'],st)+p('低R28区間中27区間が支持条件を通過。最初の0.61–0.65は固定中心であっても幅比が既存上限2を超えて未支持。高R支持点の範囲はR=2.00–4.45で、4.30–4.35に未支持の穴がある。4.45以上も未支持であり、穴をつなぐ外挿・追加veto・単一上限は設定しない。'))
bad=[r for r in summary['rows'] if not r['quality']]
pages.append(head('9. 品質、誤差と再現性')+p('支持フラグは従来の数値成立、形状誤差、identity、選択領域内残差等の条件で決まる。支持されたことはχ²/NDF≈1やモデルの物理的正しさを保証しない。未支持の全13区間もfit集に掲載する。')+table(['未支持R区間','主な原因'],[[f"{r['rLo']:.2f}–{r['rHi']:.2f}",'局所残差'+('・形状支持' if 'free_target_shape_or_identity_screen' in r['reasons'] else '')] for r in bad])+p('誤差はモデル・形状を条件とした収量誤差。集約では分母で重み付けして線形和をとる保守的な上限であり、応答形状・較正位置・背景選択・高Rで除いた他核種の完全な系統誤差ではない。背景なしでPurityが1近くなる点も仮定の影響を含む。')+p('検証：d/t/³HeのGit一致、元ヒスト不変、R<2の背景=0、高R d/t/³He成分=0を確認。556件の独立生計数、417点の³He半値包含、固定中心19件が元linear値と全桁一致し、R≥0.8の180 native区間・84集約点と全PID係数が前版と一致することを確認した。旧入力SHAと全PDFの読戻しを検証し、completed.jsonに保存する。')+p('再現ツール：tools/nuclear_pid_he4_fixedmean_20261022/。')+p('設定・係数・選択条件・集約値・各native区間の値、ソースとコマンドを同梱。Git元ヘッダと候補差分も保存。')+p('mainconf：main_auau3p85_anaFemtoLambda_4He_KFParticle.yaml。STAR SL24y/ROOT5の既存wrapperで実行。2021 P24iyキャッシュのローカルstudyで新規jobid/configlogなし。本番PID反映・job投入・hadd・commit/pushなし。'))
comparison=run.read(R/'low_comparison.json');fixed=run.read(R/'fixed_center_validation.json')
pages.append(head('10. 低R再評価の変更前後')+p('今回変えたのは0.61–0.80の19 native区間だけ。前版Gaussian＋linearの⁴He中心を各区間で固定し、linearを除去。幅と核種収量をfitし直した。PID関数の係数は変更していない。')+table(['R区間','前版linear P','中心固定・背景なし P','支持'],[[f"{r['rLo']:.2f}–{r['rHi']:.2f}",f"{100*r['previousLinearPurity']:.4f}%",f"{100*r['fixedMeanGaussianOnlyPurity']:.4f}%",'yes' if r['quality'] else 'no'] for r in comparison])+table(['R下限','固定⁴He中心','前版幅','背景なし幅'],[[f"{r['rLo']:.2f}",f"{r['mean']:.4f}",f"{r['previousLinearSigma']:.4f}",f"{r['sigma']:.4f}"] for r in fixed[:4]])+p('中心・幅の単位はkeV/cm。中心は元の値を正確に維持した一方、背景のないモデルでは⁴He幅が広がる。0.61–0.62と0.62–0.63の幅比は元のidentity参照に対して2.08694、2.12677で、既存上限2を超える。採否条件は緩めず、最初の集約点を図の下端の×として示した。')+p('linearを除いたためモデル内の混入が減り、Purityが1に近づく。これはデータ中の背景が消えたという測定ではない。中心固定のみで形状残差や幅の問題が解消したとは主張しない。')+p('固定中心のmeanError=0は条件付きfit誤差であり、同じデータから前版で推定した中心の不確かさを含まない。縦誤差も条件付き収量誤差で、較正・応答・背景モデルの完全な系統誤差ではない。'))
preamble=r'''\documentclass[10pt,a4paper]{article}
\usepackage[margin=18mm,top=20mm,bottom=18mm]{geometry}
\usepackage{fontspec,graphicx,amsmath,booktabs,xcolor,fancyhdr}
\setmainfont{DejaVu Sans}
\newfontfamily\jpfont[Path=/usr/share/fonts/google-droid-sans-fonts/,BoldFont=DroidSansFallbackFull.ttf,BoldFeatures={FakeBold=1.2}]{DroidSansFallbackFull.ttf}
\XeTeXlinebreaklocale "ja"
\XeTeXlinebreakskip=0pt plus 1pt
\definecolor{navy}{HTML}{17385D}
\pagestyle{fancy}\fancyhf{}\fancyhead[L]{\small STAR / Nuclear PID and Purity}\fancyhead[R]{\small He4 / region-specific models}\fancyfoot[C]{\thepage}\setlength{\headheight}{14pt}
\emergencystretch=2em\setlength{\parindent}{0pt}\setlength{\parskip}{5pt}\linespread{1.10}
\begin{document}
'''
tex=preamble+'\n\\clearpage\n'.join(pages)+'\n\\end{document}\n';(R/(STEM+'.tex')).write_text(tex)
for attempt in range(2):
 cmd=['xelatex','-interaction=nonstopmode','-halt-on-error','-file-line-error',STEM+'.tex']
 with (R/f'latex_build_{attempt+1}.log').open('w') as f:status=subprocess.run(cmd,cwd=R,stdout=f,stderr=subprocess.STDOUT).returncode
 assert not status,'XeLaTeX failed'
log=(R/(STEM+'.log')).read_text();assert 'Missing character:' not in log;run.save(R/'report_build.json',dict(stem=STEM,expectedPages=len(pages),missingGlyphs=False,overfullBoxes=log.count('Overfull \\hbox'),command=cmd));print('REPORT_BUILT',len(pages),'pages, overfull',log.count('Overfull \\hbox'),flush=True)
