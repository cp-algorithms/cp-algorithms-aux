#!/usr/bin/env python3
"""
Inject minified code into documentation markdown files.
This script reads minified versions from cp-algo/min/ and cp-algo/min-bundled/
and adds minifiedCode and minifiedBundledCode fields to the documentation markdown files.
"""

import os
import sys
from pathlib import Path
import re
import json
import yaml


def inject_minified_to_markdown(markdown_file, minified_code=None, minified_bundled_code=None):
    """Inject minified code into markdown file's front matter using proper YAML parsing."""
    try:
        with open(markdown_file, 'r', encoding='utf-8') as f:
            content = f.read()
        
        # Check if file has YAML front matter
        if not content.startswith('---'):
            return False
        
        # Split front matter and content at first --- and second ---
        parts = content.split('---', 2)
        if len(parts) < 3:
            return False
        
        front_matter_str = parts[1]
        body = parts[2]
        
        # Parse YAML front matter
        try:
            front_matter = yaml.safe_load(front_matter_str)
        except Exception as e:
            print(f"Error parsing YAML for {markdown_file}: {e}", file=sys.stderr)
            return False
        
        if not isinstance(front_matter, dict):
            return False
        
        updated = False
        
        # The code is nested in front_matter['data'], not at root level
        if 'data' not in front_matter or not isinstance(front_matter['data'], dict):
            return False
        
        # Add or update minifiedCode in the nested data object
        if minified_code:
            front_matter['data']['minifiedCode'] = minified_code
            updated = True
        
        # Add or update minifiedBundledCode in the nested data object
        if minified_bundled_code:
            front_matter['data']['minifiedBundledCode'] = minified_bundled_code
            updated = True
        
        if not updated:
            return False
        
        # Re-serialize YAML front matter
        # Use default_flow_style=False to keep lists as blocks, allow_unicode=True for special chars
        new_front_matter_str = yaml.dump(front_matter, default_flow_style=False, allow_unicode=True, sort_keys=False)
        
        # Write updated content
        new_content = f'---{new_front_matter_str}---{body}'
        with open(markdown_file, 'w', encoding='utf-8') as f:
            f.write(new_content)
        
        return True
    except Exception as e:
        print(f"Error processing {markdown_file}: {e}", file=sys.stderr)
        return False


def main():
    markdown_dir = Path('_jekyll')
    minified_dir = Path('cp-algo/min')
    minified_bundled_dir = Path('cp-algo/min-bundled')
    
    # If _jekyll doesn't exist, try the verify-helper path (for local testing)
    if not markdown_dir.exists():
        markdown_dir = Path('.verify-helper/markdown')
    
    if not minified_dir.exists():
        print(f"Warning: {minified_dir} does not exist", file=sys.stderr)
    
    if not minified_bundled_dir.exists():
        print(f"Warning: {minified_bundled_dir} does not exist", file=sys.stderr)
    
    if not markdown_dir.exists():
        print(f"Error: {markdown_dir} does not exist", file=sys.stderr)
        sys.exit(1)
    
    print("Injecting minified code into documentation...")
    
    count = 0
    # Find all markdown files
    for md_file in markdown_dir.rglob('*.md'):
        # Get relative path without .md extension
        rel_path = md_file.relative_to(markdown_dir)
        path_without_ext = str(rel_path)[:-3]  # Remove .md
        
        minified_code = None
        minified_bundled_code = None
        
        # Try to find corresponding minified source file
        possible_extensions = ['hpp', 'cpp', 'h']
        
        for ext in possible_extensions:
            if minified_code is None and minified_dir.exists():
                minified_file = minified_dir / f"{path_without_ext}.{ext}"
                if minified_file.exists():
                    with open(minified_file, 'r', encoding='utf-8') as f:
                        minified_code = f.read()
                    break
        
        # Try to find corresponding minified bundled file
        for ext in possible_extensions:
            if minified_bundled_code is None and minified_bundled_dir.exists():
                minified_bundled_file = minified_bundled_dir / f"{path_without_ext}.{ext}"
                if minified_bundled_file.exists():
                    with open(minified_bundled_file, 'r', encoding='utf-8') as f:
                        minified_bundled_code = f.read()
                    break
        
        # Only inject if we found at least one minified version
        if (minified_code or minified_bundled_code) and inject_minified_to_markdown(md_file, minified_code, minified_bundled_code):
            count += 1
            print(f"  Updated: {path_without_ext}")
    
    print(f"\nUpdated {count} documentation files")
    return 0


if __name__ == '__main__':
    sys.exit(main())
