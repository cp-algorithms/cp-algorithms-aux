#!/usr/bin/env python3
"""
Generate minified versions of source library files.
This script processes source files from cp-algo/ and creates
minified versions in two places:
1. .competitive-verifier/minified/ (for CI/documentation)
2. cp-algo/min/ (committed to repo for direct access)

Unlike bundled versions, minified versions preserve the original file structure
without inlining dependencies.
"""

import os
import re
import sys
import subprocess
from pathlib import Path


def minify_cpp(code):
    """Delegate to the shared minifier so CI and local minified files match."""
    minifier = Path(__file__).resolve().parents[2] / 'cp-algo/minify.py'
    try:
        result = subprocess.run(
            ['python3', str(minifier)],
            input=code,
            capture_output=True,
            text=True,
            check=True,
        )
        return result.stdout
    except FileNotFoundError:
        print(f"Warning: minifier not found at {minifier}, skipping minification", file=sys.stderr)
        return code
    except subprocess.CalledProcessError as e:
        print(f"Warning: minifier failed ({e}); skipping minification", file=sys.stderr)
        if e.stdout:
            return e.stdout
        return code


def process_file(bundled_path, minified_path, committed_path=None):
    """Process a single file and create minified version(s)."""
    try:
        with open(bundled_path, 'r', encoding='utf-8') as f:
            code = f.read()
        
        minified_code = minify_cpp(code)
        
        # Create minified version in .competitive-verifier/minified/
        minified_path.parent.mkdir(parents=True, exist_ok=True)
        with open(minified_path, 'w', encoding='utf-8') as f:
            f.write(minified_code)
        
        # Also create in committed minified/ directory if path provided
        if committed_path:
            committed_path.parent.mkdir(parents=True, exist_ok=True)
            with open(committed_path, 'w', encoding='utf-8') as f:
                f.write(minified_code)
        
        original_size = len(code)
        minified_size = len(minified_code)
        reduction = original_size - minified_size
        reduction_pct = 100 * (1 - minified_size / original_size) if original_size > 0 else 0
        
        print(f"  {bundled_path.name}: {original_size:,} → {minified_size:,} bytes (-{reduction_pct:.1f}%)")
        return True
    except Exception as e:
        print(f"  ERROR processing {bundled_path}: {e}", file=sys.stderr)
        return False


def main():
    # Source directory to minify
    source_dir = Path('cp-algo')
    bundled_dir = Path('.competitive-verifier/bundled')
    
    # Output directories
    minified_ci_dir = Path('.competitive-verifier/minified')
    minified_bundled_ci_dir = Path('.competitive-verifier/minified-bundled')
    minified_committed_dir = Path('cp-algo/min')
    
    # Verify source directory exists
    if not source_dir.exists():
        print(f"Error: {source_dir} does not exist", file=sys.stderr)
        sys.exit(1)
    
    # Clear output directories
    for d in [minified_ci_dir, minified_bundled_ci_dir]:
        if d.exists():
            import shutil
            shutil.rmtree(d)
        d.mkdir(parents=True, exist_ok=True)
    
    minified_committed_dir.mkdir(parents=True, exist_ok=True)
    
    print("Generating minified versions from source files...")
    
    total_files = 0
    processed_files = 0
    
    # Process all source files in cp-algo (but not in cp-algo/min* itself)
    for src_file in source_dir.rglob('*'):
        # Skip files in min directories
        if 'min' in src_file.parts:
            continue
        
        if src_file.is_file() and src_file.suffix in ['.hpp', '.cpp', '.h']:
            total_files += 1
            
            # Calculate relative path within cp-algo
            rel_path = src_file.relative_to(source_dir)
            
            # Output paths for minified source
            minified_ci_file = minified_ci_dir / 'cp-algo' / rel_path
            minified_committed_file = minified_committed_dir / rel_path
            
            if process_file(src_file, minified_ci_file, minified_committed_file):
                processed_files += 1
    
    print(f"\nProcessed {processed_files}/{total_files} source files")
    print(f"Generated source minified in:")
    print(f"  - .competitive-verifier/minified/cp-algo/")
    print(f"  - cp-algo/min/")
    
    # Now process bundled versions if they exist
    # NOTE: Bundled minified files are ONLY for CI/docs injection
    # They are NOT committed to the repo to avoid bloat
    if bundled_dir.exists():
        print(f"\nGenerating minified bundled versions from bundled files...")
        total_bundled = 0
        processed_bundled = 0
        
        for bundled_file in bundled_dir.rglob('*'):
            # Skip verify directory - these are test files with #line directives
            # that get corrupted during minification
            if 'verify' in bundled_file.parts:
                continue
            
            # Skip files that are already in min or min-bundled directories
            # (these would be from previous minification runs in the bundled dir)
            rel_path_str = str(bundled_file.relative_to(bundled_dir))
            # Check for 'min' or 'min-bundled' as directory components
            if any(part in ['min', 'min-bundled'] for part in bundled_file.relative_to(bundled_dir).parts):
                print(f"  Skipping already-minified: {rel_path_str}")
                continue
            
            if bundled_file.is_file() and bundled_file.suffix in ['.hpp', '.cpp', '.h']:
                total_bundled += 1
                
                # Calculate relative path within bundled
                rel_path = bundled_file.relative_to(bundled_dir)
                
                # Strip the 'cp-algo/' prefix if present (bundled dir already has full structure)
                # so we don't get cp-algo/min-bundled/cp-algo/...
                if rel_path.parts[0] == 'cp-algo':
                    rel_path = Path(*rel_path.parts[1:])
                
                # Output path for minified bundled (CI only, NOT committed)
                minified_bundled_ci_file = minified_bundled_ci_dir / rel_path
                
                # Only write to CI directory
                minified_bundled_ci_file.parent.mkdir(parents=True, exist_ok=True)
                with open(bundled_file, 'r', encoding='utf-8') as f:
                    content = f.read()
                minified = minify_cpp(content)
                with open(minified_bundled_ci_file, 'w', encoding='utf-8') as f:
                    f.write(minified)
                processed_bundled += 1
        
        print(f"\nProcessed {processed_bundled}/{total_bundled} bundled files")
        print(f"Generated bundled minified in:")
        print(f"  - .competitive-verifier/minified-bundled/ (CI only, not committed)")
    
    if processed_files < total_files and total_files > 0:
        sys.exit(1)


if __name__ == '__main__':
    main()
