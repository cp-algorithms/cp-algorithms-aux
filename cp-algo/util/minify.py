#!/usr/bin/env python3
"""
Simple C++ code minifier for competitive programming.
Removes comments, extra whitespace, and compresses the code.
"""

import re
import sys

def minify_cpp(code):
    # Remove single-line comments
    code = re.sub(r'//.*?$', '', code, flags=re.MULTILINE)
    
    # Remove multi-line comments (but preserve #line directives content)
    code = re.sub(r'/\*.*?\*/', '', code, flags=re.DOTALL)
    
    # Remove empty lines
    code = re.sub(r'\n\s*\n', '\n', code)
    
    # Remove leading/trailing whitespace on each line
    lines = [line.strip() for line in code.split('\n')]
    
    # Keep preprocessor directives on their own lines
    result = []
    for line in lines:
        if line.startswith('#'):
            result.append(line)
        elif line:
            result.append(line)
    
    code = '\n'.join(result)
    
    # Compress spaces (but not in preprocessor directives or strings)
    def compress_line(line):
        if line.startswith('#'):
            return line
        # Keep strings intact
        parts = re.split(r'("(?:[^"\\]|\\.)*")', line)
        for i in range(0, len(parts), 2):
            # Compress multiple spaces to one
            parts[i] = re.sub(r'\s+', ' ', parts[i])
            # Remove spaces around operators (carefully)
            parts[i] = re.sub(r'\s*([+\-*/%=<>!&|^~,;:?(){}[\]])\s*', r'\1', parts[i])
        return ''.join(parts)
    
    lines = [compress_line(line) for line in code.split('\n')]
    code = '\n'.join(lines)
    
    return code

if __name__ == '__main__':
    if len(sys.argv) > 1:
        with open(sys.argv[1], 'r') as f:
            code = f.read()
    else:
        code = sys.stdin.read()
    
    minified = minify_cpp(code)
    print(minified)
    
    # Print stats to stderr
    original_size = len(code)
    minified_size = len(minified)
    print(f"Original: {original_size} bytes", file=sys.stderr)
    print(f"Minified: {minified_size} bytes", file=sys.stderr)
    print(f"Reduction: {original_size - minified_size} bytes ({100*(1-minified_size/original_size):.1f}%)", file=sys.stderr)
