#!/usr/bin/env python3
"""Alternate saved oj-verify binaries on the largest cached input of each problem."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import resource
import statistics
import subprocess
import time

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--before', type=Path, required=True)
parser.add_argument('--after', type=Path, required=True)
parser.add_argument('--cache', type=Path, default=Path('.verify-helper/cache'))
parser.add_argument('--repeat', type=int, default=5)
parser.add_argument('--cases', type=int, default=1)
parser.add_argument('--output', type=Path, required=True)
parser.add_argument('tests', type=Path, nargs='+')
args = parser.parse_args()
if args.repeat < 1 or args.cases < 1:
    parser.error('--repeat and --cases must be positive')
cpu = min(os.sched_getaffinity(0))
os.sched_setaffinity(0, {cpu})
results = []
for test in args.tests:
    url = re.search(r'#define PROBLEM "([^"]+)"', test.read_text())[1]
    key = hashlib.md5(url.encode()).hexdigest()
    cases = sorted((args.cache / key / "test").glob("*.in"),
                   key=lambda p: (p.stat().st_size, p.name), reverse=True)[:args.cases]
    if not cases:
        raise FileNotFoundError(f"No cached inputs for {test}")
    for case in cases:
        tag = test.parent.name + '-' + test.stem
        binaries = [args.before / (tag + '.out'), args.after / (tag + '.out')]
        samples = [[], []]
        for rep in range(args.repeat + 1):
            for side in ([0, 1] if rep % 2 == 0 else [1, 0]):
                with case.open('rb') as src:
                    start = resource.getrusage(resource.RUSAGE_CHILDREN)
                    wall = time.monotonic()
                    subprocess.run([str(binaries[side].resolve())], stdin=src,
                                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                                   check=True)
                    wall = time.monotonic() - wall
                    end = resource.getrusage(resource.RUSAGE_CHILDREN)
                if rep:  # One untimed warm-up for each executable.
                    samples[side].append({'wall': wall, 'cpu': end.ru_utime + end.ru_stime - start.ru_utime - start.ru_stime})
        medians = [statistics.median(x['wall'] for x in side) for side in samples]
        cpu_medians = [statistics.median(x['cpu'] for x in side) for side in samples]
        row = dict(test=str(test), case=str(case), cpu=cpu, before=medians[0], after=medians[1],
                   ratio=medians[1] / medians[0], cpu_before=cpu_medians[0], cpu_after=cpu_medians[1],
                   cpu_ratio=cpu_medians[1] / cpu_medians[0], samples=samples,
                   sha256=[hashlib.sha256(p.read_bytes()).hexdigest() for p in binaries])
        results.append(row)
        args.output.write_text(json.dumps(results, indent=2) + '\n')
        print(f'{test.stem:24} {medians[0]:.4f} -> {medians[1]:.4f} s ({row["ratio"]:.3f}x wall, {row["cpu_ratio"]:.3f}x CPU)', flush=True)
