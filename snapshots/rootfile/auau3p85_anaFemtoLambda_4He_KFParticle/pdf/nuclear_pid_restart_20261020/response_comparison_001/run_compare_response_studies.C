void run_compare_response_studies() {
  if (gSystem->Load("libMinuit") < 0) { gSystem->Exit(1); return; }
  gSystem->SetBuildDir("/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_restart_20261020/response_comparison_001/build", kFALSE);
  if (!gSystem->CompileMacro("/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_restart_20261020/response_comparison_001/source/tools/nuclear_pid_restart_20261020/DrawResponseComparison.C", "k")) { gSystem->Exit(1); return; }
  Int_t error = 0;
  gROOT->ProcessLine("DrawResponseComparison(\"/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_restart_20261020/response_comparison_001/sources.tsv\",\"/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_restart_20261020/response_comparison_001/view_settings.tsv\",\"/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_restart_20261020/response_comparison_001/response_comparison.pdf\");", &error);
  if (error) gSystem->Exit(1);
}
