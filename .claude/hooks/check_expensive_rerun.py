#!/usr/bin/env python3
"""PreToolUse hook (Bash): any job submission that demultiplexes, aligns or marks duplicates goes to the user for confirmation.
(2026-09-24, user's choice: a prompt only — no registry, no name matching. Skipping work that is already done belongs in the
pipeline, zealgt docs/PLAN_pipeline.md §0 Task 2.)

It reads the Bash command plus every local script the command runs (`ssh hazel 'bash -s' < agent/suggested_script_*.sh`,
`sbatch <file>`, `bash <file>`), and the scripts those scripts run in turn. File names that are only mentioned are not followed. When that text submits work
(sbatch / nextflow run) and names an expensive step, it returns "ask": the user sees the permission prompt with the reason.
Stub runs are ignored. It never allows or denies; otherwise the normal permission flow runs."""
import json, os, re, sys

REPO = os.environ.get('CLAUDE_PROJECT_DIR') or os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
SUBMIT = re.compile(r'(?<![.\w])sbatch\b|\bnextflow\s+run\b')   # not the .sbatch file extension
STUB = re.compile(r'-stub(-run)?\b')
RUNS_JOBS = re.compile(r'(^|[;&|(]\s*)(\w+=\S*\s+)*(ssh|sbatch|nextflow|bash|sh)\b')   # git/grep/cat etc. never submit
EXPENSIVE = {
    'demultiplexing':    re.compile(r'q_nilhmm_(pool_run|bc2s3_batch2|gate2|pipeline)\.sh|--entry\s+demux\w*|\bnextflow\s+run\s+\S*main\.nf|\bcutadapt\b|\bsabre\b|NVS188B_\S+\.tar'),
    'alignment':         re.compile(r'q_nilhmm_(pool_run|bc2s3_batch2|gate2|pipeline)\.sh|\bnextflow\s+run\s+\S*main\.nf|\bminibwa\b|\bbwa(-mem2)?\s+mem\b'),
    'duplicate marking': re.compile(r'\bsamtools\s+markdup\b|\bMarkDuplicates\b'),
}

def script_paths(text):
    """Only scripts that are executed: `< file`, `sbatch [opts] file`, `bash|sh file`. A file name merely mentioned is not followed.
    Relative names are also looked up in the directories the launchers are run from (nilhmm/, nilhmm/bin/, PHG/bin/)."""
    cands = re.findall(r'<\s*([\w./~-]+\.(?:sh|sbatch))', text)
    cands += re.findall(r'(?<![.\w])(?:sbatch|bash|sh)\s+(?:--?\S+\s+)*([\w./~-]+\.(?:sh|sbatch))', text)
    seen = []
    for c in cands:
        c = os.path.expanduser(c)
        tries = [c] if os.path.isabs(c) else [os.path.join(REPO, d, c) for d in ('', 'nilhmm', 'nilhmm/bin', 'PHG/bin')]
        for p in tries:
            if os.path.isfile(p) and p not in seen: seen.append(p); break
    return seen

def gather(cmd):
    texts, done, queue = [cmd], set(), script_paths(cmd)
    for _ in range(2):                      # the command's scripts, then the launchers they call
        nxt = []
        for p in queue:
            if p in done: continue
            done.add(p)
            try: t = open(p, errors='replace').read()
            except OSError: continue
            texts.append(t); nxt += script_paths(t)
        queue = nxt
    return texts, sorted(done)

def strip_comments(t):
    return '\n'.join(l.split(' #')[0] for l in t.splitlines() if not l.lstrip().startswith('#'))

def main():
    try: cmd = json.load(sys.stdin).get('tool_input', {}).get('command', '')
    except Exception: sys.exit(0)
    if not cmd or not RUNS_JOBS.search(cmd): sys.exit(0)
    texts, files = gather(cmd)
    body = '\n'.join(strip_comments(t) for t in texts)
    if not SUBMIT.search(body) or STUB.search(body): sys.exit(0)
    steps = [s for s, rx in EXPENSIVE.items() if rx.search(body)]
    if not steps: sys.exit(0)
    where = ', '.join(os.path.relpath(p, REPO) for p in files) or 'the command itself'
    reason = (f'This submits {", ".join(steps)} ({where}). Confirm it is not a repeat of work already done '
              '(demux FASTQs in work/, CRAMs in results/*/cram, align_membench, bc2s3_realign).')
    print(json.dumps({'hookSpecificOutput': {'hookEventName': 'PreToolUse', 'permissionDecision': 'ask',
                                             'permissionDecisionReason': reason}}))

main()
