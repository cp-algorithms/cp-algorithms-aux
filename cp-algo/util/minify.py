#!/usr/bin/env python3
"""
Simple C++ code minifier for competitive programming.
Removes comments, extra whitespace, and compresses the code.
"""

import re
import sys

def minify_cpp(code):
    # Remove comments while preserving strings
    lines = code.split('\n')
    cleaned_lines = []
    
    for line in lines:
        # Remove single-line comments, but not if // is inside a string
        in_string = False
        escape_next = False
        result = []
        i = 0
        
        while i < len(line):
            if escape_next:
                result.append(line[i])
                escape_next = False
                i += 1
                continue
            
            if line[i] == '\\' and in_string:
                result.append(line[i])
                escape_next = True
                i += 1
                continue
            
            if line[i] == '"':
                result.append(line[i])
                in_string = not in_string
                i += 1
                continue
            
            if not in_string and i + 1 < len(line) and line[i:i+2] == '//':
                # Found comment outside string, discard rest of line
                break
            
            result.append(line[i])
            i += 1
        
        cleaned_lines.append(''.join(result))
    
    code = '\n'.join(cleaned_lines)
    
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
        
        # Split by string literals while keeping them
        result = []
        current = []
        i = 0
        
        while i < len(line):
            if line[i] == '"':
                # Found start of string - save accumulated code
                if current:
                    code_part = ''.join(current)
                    # Compress the code part
                    code_part = re.sub(r'\s+', ' ', code_part)
                    code_part = re.sub(r'\s*([+\-*/%=<>!&|^~,;:?(){}[\]])\s*', r'\1', code_part)
                    result.append(code_part)
                    current = []
                
                # Now collect the entire string literal
                string_chars = ['"']
                i += 1
                while i < len(line):
                    if line[i] == '\\' and i + 1 < len(line):
                        # Escape sequence
                        string_chars.append(line[i])
                        string_chars.append(line[i + 1])
                        i += 2
                    elif line[i] == '"':
                        # End of string
                        string_chars.append('"')
                        i += 1
                        break
                    else:
                        string_chars.append(line[i])
                        i += 1
                
                # Keep string literal as-is
                result.append(''.join(string_chars))
            else:
                current.append(line[i])
                i += 1
        
        # Handle any remaining code
        if current:
            code_part = ''.join(current)
            code_part = re.sub(r'\s+', ' ', code_part)
            code_part = re.sub(r'\s*([+\-*/%=<>!&|^~,;:?(){}[\]])\s*', r'\1', code_part)
            result.append(code_part)
        
        return ''.join(result)
    
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
