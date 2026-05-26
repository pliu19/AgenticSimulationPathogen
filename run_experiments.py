#!/usr/bin/python3
"""
Orchestration script: runs all 3 regimens x 36 (p, q) combinations = 108 simulations.
Uses a semaphore to cap concurrent processes at MAX_PARALLEL.
r = p, s = 0.1 * p.
"""

import glob
import os
import subprocess
import sys
import threading
import time
from itertools import product

# Paths are derived from this file's location so the script runs unchanged on
# any machine (Linux / macOS / Windows). Override with env vars if needed.
REPO_ROOT = os.path.dirname(os.path.abspath(__file__))
SRC_DIR   = os.environ.get('ABM_SRC_DIR', os.path.join(REPO_ROOT, 'src'))
MAIN      = os.environ.get('ABM_MAIN',    os.path.join(SRC_DIR, 'main.py'))


def _resolve_python():
    """Prefer the project's .venv interpreter; fall back to the current one.

    Honors ABM_PYTHON if set. Looks for a venv at <repo>/.venv (Linux/macOS
    use bin/python; Windows uses Scripts/python.exe).
    """
    override = os.environ.get('ABM_PYTHON')
    if override:
        return override
    venv = os.path.join(REPO_ROOT, '.venv')
    candidates = [
        os.path.join(venv, 'bin', 'python'),         # Unix
        os.path.join(venv, 'Scripts', 'python.exe'), # Windows
    ]
    for c in candidates:
        if os.path.isfile(c):
            return c
    print(f'[WARN] No .venv found at {venv}; using {sys.executable}. '
          f'Create one with: python -m venv .venv && '
          f'.venv/bin/pip install -r requirements.txt', flush=True)
    return sys.executable


PYTHON      = _resolve_python()
MAX_PARALLEL = 20  # ~20 concurrent processes on 28 logical cores
MONITOR_INTERVAL = 120  # seconds between monitor prints

REGIMENS  = ['fixed', 'cycling', 'mixing']
PQ_VALUES = [round(v * 0.1, 1) for v in range(1, 7)]   # 0.1 to 0.6
PQ_PAIRS  = list(product(PQ_VALUES, PQ_VALUES))          # 36 (p, q) combinations

RUN_NAME = 'run3'

SHARED_ARGS = [
    '--num_runs',    '1000',
    '--days',        '1460',
    '--empiric_hours', '72',
    '--cycle_months',  '6',
    '--run_name',    RUN_NAME,
]

total      = len(REGIMENS) * len(PQ_PAIRS)
done       = 0
done_lock  = threading.Lock()
semaphore  = threading.Semaphore(MAX_PARALLEL)
claimed    = set()   # tracks jobs already claimed to prevent double-runs
start_time = time.time()


# ─── Monitor thread ───────────────────────────────────────────────────────────

def _get_python_memory_mb():
    """Return (process_count, total_MB) across all running python processes.

    Cross-platform: uses `tasklist` on Windows and `ps` elsewhere.
    """
    try:
        if os.name == 'nt':
            out = subprocess.check_output(
                ['tasklist', '/FI', 'IMAGENAME eq python.exe', '/FO', 'CSV', '/NH'],
                text=True, stderr=subprocess.DEVNULL)
            total_kb = 0
            count = 0
            for line in out.strip().splitlines():
                parts = line.strip('"').split('","')
                if len(parts) >= 5:
                    mem_str = parts[4].replace(',', '').replace(' K', '').strip()
                    try:
                        total_kb += int(mem_str)
                        count += 1
                    except ValueError:
                        pass
            return count, total_kb / 1024
        # Linux / macOS: RSS in KB from `ps`, filter to python procs by command name.
        out = subprocess.check_output(
            ['ps', '-eo', 'comm=,rss='],
            text=True, stderr=subprocess.DEVNULL)
        total_kb = 0
        count = 0
        for line in out.strip().splitlines():
            parts = line.split(None, 1)
            if len(parts) != 2:
                continue
            comm, rss = parts[0], parts[1].strip()
            if 'python' not in os.path.basename(comm).lower():
                continue
            try:
                total_kb += int(rss)
                count += 1
            except ValueError:
                pass
        return count, total_kb / 1024
    except Exception:
        return 0, 0.0


def monitor_thread():
    """Periodically print progress, speed, ETA, and memory usage."""
    while True:
        time.sleep(MONITOR_INTERVAL)
        with done_lock:
            n_done = done
        elapsed   = time.time() - start_time
        rate      = n_done / elapsed * 3600 if elapsed > 0 else 0   # jobs/hour
        remaining = total - n_done
        eta_h     = remaining / rate if rate > 0 else float('inf')
        n_proc, mem_mb = _get_python_memory_mb()
        print(
            f'\n[MONITOR] {n_done}/{total} done | '
            f'elapsed {elapsed/3600:.1f}h | '
            f'rate {rate:.1f} jobs/h | '
            f'ETA {eta_h:.1f}h | '
            f'{n_proc} python procs | '
            f'RAM {mem_mb/1024:.1f} GB\n',
            flush=True)
        if n_done >= total:
            break


# ─── Job functions ────────────────────────────────────────────────────────────

def is_done(regimen, p_val, q_val):
    """Return True if a pkl file already exists for this run."""
    out_dir = os.path.join(SRC_DIR, 'log', RUN_NAME,
                           f'{regimen}_h72_p{p_val}_q{q_val}_r{p_val}_n64_w16')
    return bool(glob.glob(os.path.join(out_dir, '*.pkl')))


def run_one(regimen, p_val, q_val):
    global done
    r_val = p_val
    s_val = round(0.1 * p_val, 4)
    key = (regimen, p_val, q_val)
    with done_lock:
        if key in claimed or is_done(regimen, p_val, q_val):
            done += 1
            print(f'  [{done}/{total}] {regimen} p={p_val} q={q_val}: SKIPPED (already done)')
            return
        claimed.add(key)
    cmd = [
        PYTHON, MAIN,
        '--empiric_regimen', regimen,
        '--p', str(p_val),
        '--q', str(q_val),
        '--r', str(r_val),
        '--s', str(s_val),
    ] + SHARED_ARGS
    with semaphore:
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                             text=True, cwd=SRC_DIR)
        stdout, _ = p.communicate()
    with done_lock:
        done += 1
        status = 'OK' if p.returncode == 0 else f'FAILED (code {p.returncode})'
        print(f'  [{done}/{total}] {regimen} p={p_val} q={q_val}: {status}')
        if stdout.strip():
            print(f'    {stdout.strip()}')


# ─── Entry point ──────────────────────────────────────────────────────────────

monitor = threading.Thread(target=monitor_thread, daemon=True)
monitor.start()

threads = []
for p_val, q_val in PQ_PAIRS:
    for regimen in REGIMENS:
        t = threading.Thread(target=run_one, args=(regimen, p_val, q_val))
        t.start()
        threads.append(t)

for t in threads:
    t.join()

print('\nAll runs complete.')
