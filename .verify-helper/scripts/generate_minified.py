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
from pathlib import Path


def minify_cpp(code):
    """Minify C++ code while preserving header guards and removing unnecessary whitespace."""
    lines = code.split('\n')
    result = []
    in_multiline_comment = False
    header_guard = None
    endif_lines = []
    
    for i, line in enumerate(lines):
        stripped = line.strip()
        
        # Handle multiline comments
        if '/*' in stripped:
            in_multiline_comment = True
        if '*/' in stripped:
            in_multiline_comment = False
            continue
        if in_multiline_comment:
            continue
        
        # Remove single-line comments
        if '//' in stripped:
            stripped = stripped.split('//')[0].strip()
        
        # Skip empty lines
        if not stripped:
            continue
        
        # Detect and preserve header guards
        if stripped.startswith('#ifndef '):
            header_guard = stripped
            result.append(stripped)
            continue
        elif stripped.startswith('#define ') and header_guard and len(result) == 1:
            result.append(stripped)
            continue
        elif stripped == '#endif' and header_guard:
            endif_lines.append(stripped)
            continue
        
        # Keep preprocessor directives as-is (they need to be on own lines)
        if stripped.startswith('#'):
            result.append(stripped)
            continue
        
        # Compress spaces in code lines (but preserve strings)
        def compress_line(line):
            # Split by strings to preserve their content
            parts = re.split(r'("(?:[^"\\]|\\.)*")', line)
            for i in range(0, len(parts), 2):
                # Compress multiple spaces to one
                parts[i] = re.sub(r'\s+', ' ', parts[i])
                # Remove spaces around operators (but keep space for clarity in some cases)
                parts[i] = re.sub(r'\s*([+\-*/%=<>!&|^~,;:?(){}[\]])\s*', r'\1', parts[i])
                # Remove leading/trailing spaces
                parts[i] = parts[i].strip()
            return ''.join(parts)
        
        compressed = compress_line(stripped)
        if compressed:
            result.append(compressed)
    
    # Join lines - try to put on same line when possible
    code = '\n'.join(result)
    
    # Remove newlines after opening braces and before closing braces
    code = re.sub(r'\{\n', '{', code)
    code = re.sub(r'\n\}', '}', code)
    
    # Remove newlines around colons in class/struct definitions
    code = re.sub(r'\n:', ':', code)
    code = re.sub(r':\n', ':', code)
    
    # Remove multiple consecutive newlines (keep max 1)
    code = re.sub(r'\n\n+', '\n', code)
    
    # Add back endif if we had header guards
    if endif_lines:
        code = code + '\n' + '\n'.join(endif_lines)
    
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
    minified_bundled_committed_dir = Path('cp-algo/min-bundled')
    
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
    
    for d in [minified_committed_dir, minified_bundled_committed_dir]:
        d.mkdir(parents=True, exist_ok=True)
    
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
            if '/min/' in rel_path_str or '/min-bundled/' in rel_path_str:
                continue
            
            if bundled_file.is_file() and bundled_file.suffix in ['.hpp', '.cpp', '.h']:
                total_bundled += 1
                
                # Calculate relative path within bundled
                rel_path = bundled_file.relative_to(bundled_dir)
                
                # Output paths for minified bundled
                minified_bundled_ci_file = minified_bundled_ci_dir / rel_path
                minified_bundled_committed_file = minified_bundled_committed_dir / rel_path
                
                if process_file(bundled_file, minified_bundled_ci_file, minified_bundled_committed_file):
                    processed_bundled += 1
        
        print(f"\nProcessed {processed_bundled}/{total_bundled} bundled files")
        print(f"Generated bundled minified in:")
        print(f"  - .competitive-verifier/minified-bundled/")
        print(f"  - cp-algo/min-bundled/")
    
    if processed_files < total_files and total_files > 0:
        sys.exit(1)


if __name__ == '__main__':
    main()
