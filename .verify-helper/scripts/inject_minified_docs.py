#!/usr/bin/env python3
"""
Inject minified code into documentation markdown files.
This script reads minified versions from .competitive-verifier/minified/ and
adds minifiedCode fields to the documentation markdown files.
"""

import os
import sys
from pathlib import Path
import re
import json


def get_minified_code(file_path, minified_dir):
    """Get minified code for a given file path."""
    minified_file = minified_dir / file_path
    if minified_file.exists():
        try:
            with open(minified_file, 'r', encoding='utf-8') as f:
                code = f.read()
            # Escape for YAML
            code = code.replace('\\', '\\\\')
            code = code.replace('"', '\\"')
            code = code.replace('\n', '\\n')
            return code
        except Exception as e:
            print(f"Error reading {minified_file}: {e}", file=sys.stderr)
    return None


def inject_minified_to_markdown(markdown_file, minified_code):
    """Inject minified code into markdown file's front matter."""
    try:
        with open(markdown_file, 'r', encoding='utf-8') as f:
            content = f.read()
        
        # Check if file has YAML front matter
        if not content.startswith('---'):
            return False
        
        # Split front matter and content
        parts = content.split('---', 2)
        if len(parts) < 3:
            return False
        
        front_matter = parts[1]
        body = parts[2]
        
        # Parse front matter to check if minifiedCode already exists
        if 'minifiedCode:' in front_matter:
            # Replace existing minifiedCode
            front_matter = re.sub(
                r'  minifiedCode: ".*?"(?=\n  [a-zA-Z_]|\n$)',
                f'  minifiedCode: "{minified_code}"',
                front_matter,
                flags=re.DOTALL
            )
        else:
            # Add minifiedCode after bundledCode if it exists
            if 'bundledCode:' in front_matter:
                front_matter = re.sub(
                    r'(  bundledCode: ".*?")(\n  [a-zA-Z_]|\n$)',
                    rf'\1\n  minifiedCode: "{minified_code}"\2',
                    front_matter,
                    flags=re.DOTALL
                )
            else:
                # Add at the end of front matter
                front_matter = front_matter.rstrip() + f'\n  minifiedCode: "{minified_code}"'
        
        # Write updated content
        new_content = f'---{front_matter}---{body}'
        with open(markdown_file, 'w', encoding='utf-8') as f:
            f.write(new_content)
        
        return True
    except Exception as e:
        print(f"Error processing {markdown_file}: {e}", file=sys.stderr)
        return False


def main():
    markdown_dir = Path('.verify-helper/markdown')
    minified_dir = Path('cp-algo/min')
    
    if not minified_dir.exists():
        print(f"Warning: {minified_dir} does not exist", file=sys.stderr)
        return 0
    
    if not markdown_dir.exists():
        print(f"Error: {markdown_dir} does not exist", file=sys.stderr)
        sys.exit(1)
    
    print("Injecting minified code into documentation...")
    
    count = 0
    # Find all markdown files with bundledCode
    for md_file in markdown_dir.rglob('*.md'):
        # Get relative path without .md extension
        rel_path = md_file.relative_to(markdown_dir)
        path_without_ext = str(rel_path)[:-3]  # Remove .md
        
        # Try to find corresponding minified file
        possible_extensions = ['hpp', 'cpp', 'h']
        minified_code = None
        
        for ext in possible_extensions:
            minified_file = minified_dir / f"{path_without_ext}.{ext}"
            if minified_file.exists():
                with open(minified_file, 'r', encoding='utf-8') as f:
                    code = f.read()
                # Escape for YAML
                code = code.replace('\\', '\\\\')
                code = code.replace('"', '\\"')
                code = code.replace('\n', '\\n')
                minified_code = code
                break
        
        if minified_code and inject_minified_to_markdown(md_file, minified_code):
            count += 1
            print(f"  Updated: {path_without_ext}")
    
    print(f"\nUpdated {count} documentation files")
    return 0


if __name__ == '__main__':
    sys.exit(main())
