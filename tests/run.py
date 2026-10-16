#!/usr/bin/env python3
"""Compile and run local C++ regressions, optionally against freshly minified headers."""

import argparse
from concurrent.futures import ThreadPoolExecutor
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import tempfile


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("tests", nargs="*", type=Path)
    parser.add_argument("--compiler", default=os.environ.get("CXX", "g++"))
    parser.add_argument("--jobs", type=int, default=2)
    parser.add_argument("--minified", action="store_true")
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    sources = sorted(args.tests or (root / "tests").glob("*.cpp"))
    if not sources or args.jobs < 1:
        parser.error("select at least one test and a positive job count")
    print(f"Running {len(sources)} local C++ tests ({'minified' if args.minified else 'source'} headers)", flush=True)
    with tempfile.TemporaryDirectory(prefix="cp-algo-tests-") as directory:
        build = Path(directory)
        if args.minified:
            # Generate only in the disposable tree; never modify committed cp-algo/min.
            shutil.copytree(root / "cp-algo", build / "cp-algo", ignore=shutil.ignore_patterns("min", "min-bundled"))
            for header in list((build / "cp-algo").rglob("*.hpp")):
                output = build / "cp-algo/min" / header.relative_to(build / "cp-algo")
                output.parent.mkdir(parents=True, exist_ok=True)
                result = subprocess.run(
                    ["python3", str(root / "cp-algo/minify.py")],
                    input=header.read_text(), capture_output=True, text=True, check=True,
                )
                output.write_text(result.stdout)

        def run(item):
            index, source = item
            try:
                original = source.resolve()
                if args.minified:
                    source = build / f"test-{index}.cpp"
                    source.write_text(re.sub(r'"cp-algo/(?!min/)', '"cp-algo/min/', original.read_text()))
                binary = build / f"test-{index}"
                command = [args.compiler, "-std=c++23", "-O2", *shlex.split(os.environ.get("CXXFLAGS", "")),
                           "-UNDEBUG", "-I" + str(build), "-I" + str(root), "-I" + str(original.parent)]
                if args.minified:
                    command.append("-DTEST_MINIFIED")
                result = subprocess.run(command + [str(source), "-o", str(binary)],
                                        capture_output=True, text=True, timeout=300)
                if result.returncode == 0:
                    result = subprocess.run([str(binary)], capture_output=True, text=True, timeout=120)
                return original.name, result.returncode == 0, result.stdout + result.stderr
            except (OSError, subprocess.TimeoutExpired) as error:
                return source.name, False, str(error)

        with ThreadPoolExecutor(max_workers=args.jobs) as pool:
            results = list(pool.map(run, enumerate(sources)))
        for name, passed, output in results:
            print(f"{'PASS' if passed else 'FAIL'} {name}")
            if not passed:
                print(output)
        passed = sum(result[1] for result in results)
        print(f"Passed {passed}/{len(results)} local C++ tests")
        return int(passed != len(results))


if __name__ == "__main__":
    raise SystemExit(main())
