#!/usr/bin/env python3
"""Wait for one SUMS run, validate/hadd its exact inputs, then optionally email.

Start in tmux on the submit host. This never resubmits, modifies runmeta, skips
corrupt ROOT files, overwrites an existing merge, or retries an uncertain merge.
"""
import argparse
import collections
import datetime
import email.message
import email.policy
import email.utils
import fcntl
import getpass
import hashlib
import json
import math
import os
from pathlib import Path
import re
import socket
import subprocess
import sys
import time
import xml.etree.ElementTree as ET

# Also works when imported by a test outside script/.
sys.path.insert(0, '/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/script')
import watch_job_notify as notify
import archive_partial as archive_helper


# ROOT 5 / CINT-compatible. Read every keyed object and every TTree entry. Store
# only histogram fingerprints, not histograms, to bound validation memory.
ROOT_AUDIT = r'''#include <TFile.h>
#include <TKey.h>
#include <TDirectory.h>
#include <TH1.h>
#include <TTree.h>
#include <TError.h>
#include <TSystem.h>
#include <TMath.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <string>
#include <cmath>
#include <cstdio>

bool mergeAuditError = false;
void mergeAuditHandler(Int_t level, Bool_t abort, const char* location, const char* msg) {
  if (level >= kError) mergeAuditError = true;
  DefaultErrorHandler(level, abort, location, msg);
}
struct MergeAuditHist {
  std::string schema;
  double entries, sum, weighted;
  MergeAuditHist(): entries(0), sum(0), weighted(0) {}
};
bool mergeAuditDir(TDirectory* dir, const std::string& prefix,
                   std::map<std::string, MergeAuditHist>& totals, int& objects) {
  TIter next(dir->GetListOfKeys());
  TKey* key = 0;
  while ((key = (TKey*)next())) {
    TObject* obj = key->ReadObj();
    if (!obj) return false;
    ++objects;
    std::string name = prefix + key->GetName() + ";";
    char cycle[32]; sprintf(cycle, "%d", key->GetCycle()); name += cycle;
    if (name.find('\t') != std::string::npos || name.find('\n') != std::string::npos) return false;
    if (obj->InheritsFrom(TDirectory::Class())) {
      // Directories belong to the TFile; Close() releases them.
      if (!mergeAuditDir((TDirectory*)obj, name + "/", totals, objects)) return false;
    } else {
      if (obj->InheritsFrom(TH1::Class())) {
        TH1* h = (TH1*)obj;
        char shape[128]; sprintf(shape, ":%d:%d:%d:%d", h->GetDimension(),
            h->GetNbinsX(), h->GetNbinsY(), h->GetNbinsZ());
        std::string schema = std::string(h->ClassName()) + shape;
        MergeAuditHist& sum = totals[name];
        if (!sum.schema.empty() && sum.schema != schema) { delete obj; return false; }
        sum.schema = schema;
        sum.entries += h->GetEntries();
        const bool profile = h->InheritsFrom("TProfile") || h->InheritsFrom("TProfile2D")
                             || h->InheritsFrom("TProfile3D");
        int cells = h->GetNbinsX() + 2;
        if (h->GetDimension() > 1) cells *= h->GetNbinsY() + 2;
        if (h->GetDimension() > 2) cells *= h->GetNbinsZ() + 2;
        for (int bin = 0; bin < cells; ++bin) {
          double v = h->GetBinContent(bin);
          if (!TMath::Finite(v)) { delete obj; return false; }
          // Profile bin means are not additive; compare their entries only.
          if (!profile) { sum.sum += v; sum.weighted += v * (bin + 1.0); }
        }
        if (!TMath::Finite(sum.entries) || !TMath::Finite(sum.sum)
            || !TMath::Finite(sum.weighted)) { delete obj; return false; }
      }
      if (obj->InheritsFrom(TTree::Class())) {
        TTree* t = (TTree*)obj;
        for (Long64_t e = 0; e < t->GetEntries(); ++e)
          if (t->GetEntry(e) < 0) { delete obj; return false; }
      }
      delete obj;
    }
    if (mergeAuditError) return false;
  }
  return true;
}
void merge_notify_audit(const char* listPath, const char* reportPath) {
  TH1::AddDirectory(kFALSE);
  SetErrorHandler(mergeAuditHandler);
  std::ifstream list(listPath);
  std::ofstream report(reportPath);
  if (!list.good() || !report.good()) { gSystem->Exit(61); return; }
  std::map<std::string, MergeAuditHist> totals;
  std::string path;
  int files = 0;
  while (std::getline(list, path)) {
    if (path.empty()) { gSystem->Exit(62); return; }
    std::cout << "AUDIT " << path << std::endl;
    TFile* f = TFile::Open(path.c_str(), "READ");
    if (!f || f->IsZombie() || f->TestBit(TFile::kRecovered) || mergeAuditError) {
      std::cerr << "AUDIT FAILED " << path << std::endl; gSystem->Exit(63); return;
    }
    int objects = 0;
    bool ok = mergeAuditDir(f, "", totals, objects);
    f->Close(); delete f;
    if (!ok || !objects || mergeAuditError) {
      std::cerr << "AUDIT FAILED " << path << std::endl; gSystem->Exit(64); return;
    }
    report << "FILE\t" << path << "\t" << objects << "\n";
    ++files;
  }
  if (!files || totals.empty()) { gSystem->Exit(65); return; }
  report << std::setprecision(17);
  for (std::map<std::string, MergeAuditHist>::const_iterator i = totals.begin(); i != totals.end(); ++i)
    report << "HIST\t" << i->first << "\t" << i->second.schema << "\t"
           << i->second.entries << "\t" << i->second.sum << "\t" << i->second.weighted << "\n";
  report << "DONE\t" << files << "\n";
  report.close();
  if (!report.good()) { gSystem->Exit(66); return; }
  gSystem->Exit(0);
}
'''


def root_layout(xml, jobid):
    output = xml.find("output").attrib
    if not output["toURL"].startswith("file:"):
        raise ValueError("Only a local file: destination is supported")
    directory = Path(output["toURL"][5:])
    template = Path(output["fromScratch"]).name
    if (not directory.is_absolute() or template.count("$JOBID") != 1
            or not template.endswith("_$JOBID.root")):
        raise ValueError("Expected an absolute directory and *_$JOBID.root output template")
    # /star/u aliases may exist only on the host; the batch container binds the
    # canonical GPFS tree. Use that path for inputs, staging, and final output.
    directory = directory.resolve()
    base = template.replace("$JOBID", jobid)
    return directory, template, directory / (Path(base).stem + "_merge.root")


def select_inputs(xml, jobid, indices, output, outcomes, allow_partial):
    directory, template, destination = root_layout(xml, jobid)
    excluded = set(output["missingRoots"] + output["badRootHeaders"]
                   + output["analysisCompletionUnconfirmed"])
    partial = bool(excluded or outcomes != {"exit_0": len(indices)})
    # A malformed ROOT file is never silently skipped, even in partial mode.
    if output["badRootHeaders"]:
        raise ValueError("ROOT headers are corrupt: {}".format(output["badRootHeaders"]))
    if partial and not allow_partial:
        raise ValueError("Some jobs/outputs are incomplete; --allow-partial was not given")
    selected = [i for i in indices if i not in excluded]
    if not selected:
        raise ValueError("No completed, valid-header subjob ROOT files to merge")
    paths = [directory / template.replace("$JOBID", "{}_{}".format(jobid, i)) for i in selected]
    return paths, selected, sorted(excluded), partial, destination


def fingerprints(paths):
    result = []
    for path in paths:
        if any(char in str(path) for char in ("\n", "\r", "\t")):
            raise ValueError("Unsupported control character in ROOT path")
        stat = path.stat()
        result.append(dict(path=str(path), size=stat.st_size, mtimeNs=stat.st_mtime_ns,
                           inode=stat.st_ino, device=stat.st_dev))
    return result


def run_logged(command, log):
    with log.open("x") as out:
        out.write("COMMAND " + json.dumps(command) + "\n")
        out.flush()
        result = subprocess.run(command, stdout=out, stderr=subprocess.STDOUT)
    if result.returncode:
        raise RuntimeError("Command failed ({}) — {}".format(result.returncode, log))
    # ROOT 5 can emit read/write errors without returning a non-zero status.
    with log.open(errors="replace") as source:
        if any(re.search(r"(?:Error|Fatal) in <|\bfatal error:|\bsegmentation (?:violation|fault)",
                         line, re.IGNORECASE) for line in source):
            raise RuntimeError("ROOT error in command log: {}".format(log))


def write_list(path, paths):
    with path.open("x") as out:
        out.write("".join(str(p) + "\n" for p in paths))
        out.flush()
        os.fsync(out.fileno())


def audit(paths, prefix, state_dir, runner):
    list_path = state_dir / (prefix + ".list")
    report = state_dir / (prefix + ".tsv")
    write_list(list_path, paths)
    macro = state_dir / "merge_notify_audit.C"
    # Compile the custom std::map value type; ROOT 5 CINT has no dictionary for it.
    invocation = "{}+({},{})".format(macro, json.dumps(str(list_path)), json.dumps(str(report)))
    run_logged(runner + ["root4star", "-b", "-q", invocation], state_dir / (prefix + ".log"))
    rows = report.read_text().splitlines()
    if not rows or rows[-1] != "DONE\t{}".format(len(paths)):
        raise RuntimeError("ROOT audit report is incomplete: {}".format(report))
    found_files = [r.split("\t")[1] for r in rows if r.startswith("FILE\t")]
    if found_files != [str(p) for p in paths]:
        raise RuntimeError("ROOT audit input list mismatch")
    histograms = {}
    for row in rows:
        fields = row.split("\t")
        if fields[0] == "HIST":
            if len(fields) != 6 or fields[1] in histograms:
                raise RuntimeError("Invalid histogram audit record")
            values = tuple(float(v) for v in fields[3:])
            if not all(math.isfinite(v) for v in values):
                raise RuntimeError("Non-finite histogram summary")
            histograms[fields[1]] = (fields[2], values)
    if not histograms:
        raise RuntimeError("No histograms were validated")
    return histograms


def compare_histograms(expected, actual):
    if set(expected) != set(actual):
        raise RuntimeError("Merged histogram names do not match input union")
    for name, (schema, values) in expected.items():
        other_schema, other_values = actual[name]
        if schema != other_schema or not all(math.isclose(a, b, rel_tol=1e-7, abs_tol=1e-7)
                                            for a, b in zip(values, other_values)):
            raise RuntimeError("Merged histogram is not the input sum: {}".format(name))


def publish_no_clobber(source, destination):
    # Staging resides on the same filesystem. link() fails atomically if an
    # output already exists, unlike rename()/replace(). Retain staging for audit.
    os.link(str(source), str(destination))


def merge_once(paths, destination, state_dir, runner, batch_size):
    attempt = state_dir / "merge_attempt.json"
    if attempt.exists():
        raise RuntimeError("Merge was already attempted; inspect records before any manual recovery")
    if os.path.lexists(str(destination)):
        raise FileExistsError("Refusing to overwrite existing merge: {}".format(destination))
    before = fingerprints(paths)
    notify.save_json(attempt, dict(startedAt=notify.now(), destination=str(destination),
                                  inputs=before, runner=runner, batchSize=batch_size))
    staging = destination.parent / (".merge_staging_" + state_dir.name)
    staging.mkdir()  # Refuse an old/inconsistent staging tree.
    notify.save_json(state_dir / "input_manifest.json", dict(inputs=before))
    write_list(state_dir / "inputs.list", paths)
    (state_dir / "merge_notify_audit.C").write_text(ROOT_AUDIT)
    expected = audit(paths, "input_audit", state_dir, runner)
    current = list(paths)
    round_number = 0
    while True:
        next_round = []
        for begin in range(0, len(current), batch_size):
            target = staging / "round_{}_part_{}.root".format(round_number, begin // batch_size)
            log = state_dir / "hadd_round_{}_part_{}.log".format(round_number, begin // batch_size)
            # No -f (overwrite) or -k (skip corrupt inputs).
            run_logged(runner + ["hadd", str(target)]
                       + [str(p) for p in current[begin:begin + batch_size]], log)
            if not target.is_file() or target.stat().st_size < 4:
                raise RuntimeError("hadd produced no ROOT output: {}".format(target))
            next_round.append(target)
        if len(next_round) == 1:
            merged = next_round[0]
            break
        current = next_round
        round_number += 1
    actual = audit([merged], "merged_audit", state_dir, runner)
    compare_histograms(expected, actual)
    if fingerprints(paths) != before:
        raise RuntimeError("An input ROOT changed during validation/merge; output not published")
    publish_no_clobber(merged, destination)
    return dict(path=str(destination), size=destination.stat().st_size,
                fingerprint=archive_helper.regular(destination),
                histogramCount=len(actual), validation="ROOT recursive read + additive histogram fingerprints",
                inputCount=len(paths), publishedAt=notify.now(), stagingDirectory=str(staging))


def make_message(recipient, meta, result):
    success = result["phase"] == "merged"
    archive = result.get("archive", {})
    label = ("hadd完了" if success and not result.get("partial")
             and archive.get("phase") != "failed" else "hadd終了・要確認")
    msg = email.message.EmailMessage(policy=email.policy.SMTP)
    msg["From"] = "{}@{}".format(getpass.getuser(), socket.getfqdn())
    msg["To"] = recipient
    msg["Subject"] = "[STAR {}] {}".format(label, meta["anaName"])
    msg["Date"] = email.utils.format_datetime(datetime.datetime.now().astimezone())
    msg["Message-ID"] = "<star-{}-merge-completion@{}>".format(meta["jobid"], socket.getfqdn())
    msg.set_content(
        "STAR解析の全job終了後のhadd処理結果です。\n\n"
        "解析: {ana}\nSUMS ID: {jobid}\n時刻（日本時間）: {when}\n"
        "状態: {phase}\nCondor終了内訳: {outcomes}\n"
        "統合入力: {selected}/{expected}\n部分統合: {partial}\n"
        "除外したsubjob番号: {excluded}\n"
        "統合ROOT: {output}\n検査: {validation}\n"
        "subjob ROOT整理: {archive_phase}\n保存先: {archive_dir}\n整理エラー: {archive_error}\n"
        "エラー: {error}\n監視・統合記録: {state}\n\n"
        "元のROOTファイルは保持しています。自動再投入は行っていません。\n"
        "ROOT読み取り検査とヒストグラム加算検査は物理QAの代わりではありません。\n".format(
            ana=meta["anaName"], jobid=meta["jobid"], when=result["checkedAt"],
            phase=result["phase"], outcomes=result.get("outcomes", {}),
            selected=len(result.get("selectedIndices", [])), expected=result["expectedSubjobs"],
            partial=result.get("partial", "not merged"), excluded=result.get("excludedIndices", []),
            output=result.get("merge", {}).get("path", "未作成"),
            validation=result.get("merge", {}).get("validation", "未完了"),
            archive_phase=archive.get("phase", "未実施"),
            archive_dir=archive.get("directory", "元の出力先に保持"),
            archive_error=archive.get("error", "なし"),
            error=result.get("error", "なし"), state=result["stateDirectory"]))
    return msg


def finish(result, meta, args, state_dir):
    if args.email:
        msg = make_message(args.email, meta, result)
        (state_dir / "message.eml").write_bytes(msg.as_bytes())
        notify.send_once(state_dir, msg, result, args.email)
        print("{} Mail accepted by local sendmail".format(notify.now()), flush=True)
    return 0 if (result["phase"] == "merged"
                 and result.get("archive", {}).get("phase") != "failed") else 1


def archive_after_merge(result, runmeta, state_dir):
    """Archive only after a durable validated merge result; never relabel hadd."""
    if result.get("phase") != "merged":
        return
    try:
        result["archive"] = archive_helper.archive_merged_subjobs(runmeta)
    except Exception as error:
        result["archive"] = dict(phase="failed", error=str(error), checkedAt=notify.now())
        print("Merged ROOT retained; subjob archive needs attention: {}".format(error), flush=True)
    notify.save_json(state_dir / "result.json", result)
    notify.save_json(state_dir / "status.json", result)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runmeta", type=Path, required=True)
    parser.add_argument("--event-log", type=Path, required=True)
    parser.add_argument("--root-runner", type=Path, required=True,
                        help="Bash wrapper invoked as runner mainconf {hadd|root4star} arguments...")
    parser.add_argument("--mainconf", type=Path)
    mail = parser.add_mutually_exclusive_group(required=True)
    mail.add_argument("--email")
    mail.add_argument("--no-email", action="store_true")
    parser.add_argument("--allow-partial", action="store_true",
                        help="After every job terminates, merge only outputs with successful analysis logs")
    parser.add_argument("--poll-sec", type=int, default=300)
    parser.add_argument("--settle-sec", type=int, default=300)
    parser.add_argument("--batch-size", type=int, default=100)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--check", action="store_true", help="Read-only status; never merge or send")
    mode.add_argument("--watch", action="store_true")
    args = parser.parse_args()
    if args.email and not re.fullmatch(r"[A-Za-z0-9.!#$%&'*+/=?^_`{|}~-]+@[A-Za-z0-9.-]+\.[A-Za-z]{2,}", args.email):
        parser.error("Invalid email address")
    if args.poll_sec < 30 or args.settle_sec < 0 or not 2 <= args.batch_size <= 1000:
        parser.error("poll-sec>=30, settle-sec>=0, and batch-size in [2,1000] are required")
    meta = json.loads(args.runmeta.read_text())
    if not re.fullmatch(r"[0-9a-fA-F]{32}", meta["jobid"]):
        parser.error("Invalid SUMS id")
    if socket.getfqdn() != meta["submit"]["host"]:
        parser.error("Run on the original submit host")
    root = Path(meta["projectRoot"])
    pattern = re.compile(re.escape(meta["jobPrefix"]) + r"_(\d+)\.list$")
    indices = sorted(int(pattern.fullmatch(p.name).group(1))
                     for p in (root / "job/run").iterdir() if pattern.fullmatch(p.name))
    if not indices or indices != list(range(len(indices))):
        parser.error("Missing/non-contiguous SUMS subjob lists")
    xml = ET.parse(meta["submit"]["joblistSnapshot"]["path"]).getroot()
    mainconf = args.mainconf or Path(meta["submit"]["mainconf"])
    if not mainconf.is_absolute():
        mainconf = root / mainconf
    mainconf = mainconf.resolve()
    runner = ["bash", str(args.root_runner.resolve()), str(mainconf)]
    if not args.root_runner.is_file() or not mainconf.is_file():
        parser.error("root-runner/mainconf not found")
    state_dir = Path(__file__).resolve().parent
    if args.watch:
        state_dir.mkdir(parents=True, exist_ok=True)
        lock = (state_dir / "monitor.lock").open("a")
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        if (state_dir / "sent.json").exists():
            print("Already notified; no duplicate merge/email", flush=True)
            return 0
        if (state_dir / "delivery_attempt.json").exists():
            raise SystemExit("Delivery outcome uncertain; inspect mail records before restarting")
        if (state_dir / "result.json").exists():
            result = json.loads((state_dir / "result.json").read_text())
            archive_after_merge(result, args.runmeta, state_dir)
            return finish(result, meta, args, state_dir)
        if (state_dir / "merge_attempt.json").exists():
            raise SystemExit("Merge outcome uncertain; inspect staging/records before restarting")
        notify.save_json(state_dir / "settings.json", dict(
            runmeta=str(args.runmeta.resolve()), eventLog=str(args.event_log.resolve()),
            email=args.email, allowPartial=args.allow_partial, runner=runner,
            pollSeconds=args.poll_sec, settleSeconds=args.settle_sec,
            batchSize=args.batch_size, expectedSubjobs=len(indices),
            pid=os.getpid(), host=socket.getfqdn(), startedAt=notify.now()))
    settled_since = None
    while True:
        try:
            submitted, terminal = notify.parse_events(args.event_log.read_text())
            queue = notify.queue_states(meta["jobid"], meta["submit"]["host"])
            complete = notify.ready(len(indices), submitted, terminal, queue)
            status = dict(checkedAt=notify.now(), expectedSubjobs=len(indices),
                          submitted=len(submitted), terminal=len(terminal), queued=len(queue),
                          queueStates=dict(collections.Counter(queue.values())),
                          outcomes=dict(collections.Counter(terminal.values())),
                          stateDirectory=str(state_dir), allTerminal=complete, phase="waiting")
            if args.check:
                print(json.dumps(status, ensure_ascii=False, indent=2))
                return 0
            if complete:
                if settled_since is None:
                    settled_since = time.monotonic()
                status["phase"] = "settling"
                if time.monotonic() - settled_since >= args.settle_sec:
                    break
            else:
                settled_since = None
            notify.save_json(state_dir / "status.json", status)
            print(json.dumps(status, ensure_ascii=False), flush=True)
        except (OSError, ValueError, RuntimeError, subprocess.SubprocessError) as error:
            if args.check:
                raise
            settled_since = None
            notify.save_json(state_dir / "status.json", dict(
                checkedAt=notify.now(), phase="waiting_retry", error=str(error),
                stateDirectory=str(state_dir)))
            print("{} Retry without merging/sending: {}".format(notify.now(), error), flush=True)
        time.sleep(args.poll_sec)
    try:
        status["phase"] = "validating_inputs"
        notify.save_json(state_dir / "status.json", status)
        output = notify.output_summary(xml, meta["jobid"], indices)
        status["output"] = output
        # Keep missing/unconfirmed indices visible even if strict selection
        # rejects this run before a merge input list can be constructed.
        status["excludedIndices"] = sorted(set(
            output["missingRoots"] + output["badRootHeaders"]
            + output["analysisCompletionUnconfirmed"]))
        paths, selected, excluded, partial, destination = select_inputs(
            xml, meta["jobid"], indices, output, status["outcomes"], args.allow_partial)
        status.update(selectedIndices=selected, excludedIndices=excluded, partial=partial,
                      phase="merging", plannedOutput=str(destination))
        notify.save_json(state_dir / "status.json", status)
        status["merge"] = merge_once(paths, destination, state_dir, runner, args.batch_size)
        status["phase"] = "merged"
    except Exception as error:
        status.update(phase="failed", error=str(error))
        print("Merge stopped: {}".format(error), flush=True)
    status["checkedAt"] = notify.now()
    notify.save_json(state_dir / "result.json", status)
    notify.save_json(state_dir / "status.json", status)
    archive_after_merge(status, args.runmeta, state_dir)
    print(json.dumps(status, ensure_ascii=False), flush=True)
    return finish(status, meta, args, state_dir)


if __name__ == "__main__":
    sys.exit(main())
