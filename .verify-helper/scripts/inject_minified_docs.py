#!/usr/bin/env python3
"""
Inject minified code into documentation markdown files.
This script reads minified versions from cp-algo/min/ and cp-algo/min-bundled/
and adds minifiedCode and minifiedBundledCode fields to the documentation markdown files.
For test files, generates minified versions on-the-fly without storing them.
"""

import os
import sys
from pathlib import Path
import re
import json
import yaml
import subprocess


class _LiteralDumper(yaml.SafeDumper):
    """Dump multiline strings using literal style so code stays readable."""


def _str_presenter(dumper, data):
    if isinstance(data, str) and "\n" in data:
        return dumper.represent_scalar("tag:yaml.org,2002:str", data, style='|')
    return dumper.represent_scalar("tag:yaml.org,2002:str", data)


# Register custom string representer
_LiteralDumper.add_representer(str, _str_presenter)


def minify_code(code):
    """Minify C++ code on-the-fly using the minify utility."""
    try:
        result = subprocess.run(
            ['python3', 'cp-algo/util/minify.py'],
            input=code,
            capture_output=True,
            text=True,
            check=True
        )
        return result.stdout
    except subprocess.CalledProcessError as e:
        print(f"Warning: Failed to minify code: {e}", file=sys.stderr)
        return None
    except FileNotFoundError:
        print("Warning: minify.py not found, skipping minification", file=sys.stderr)
        return None


def bundle_and_minify(source_file):
    """Bundle and minify a test file using oj-bundle and minify.py."""
    try:
        # First bundle with oj-bundle
        bundle_result = subprocess.run(
            ['oj-bundle', str(source_file)],
            capture_output=True,
            text=True,
            check=True
        )
        bundled_code = bundle_result.stdout
        
        # Then minify
        minified = minify_code(bundled_code)
        return minified
    except subprocess.CalledProcessError as e:
        print(f"Warning: Failed to bundle/minify {source_file}: {e}", file=sys.stderr)
        return None
    except FileNotFoundError:
        print("Warning: oj-bundle not found, skipping bundled minification", file=sys.stderr)
        return None


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
        
        # The code is nested in front_matter['data'], but create 'data' if it doesn't exist
        if 'data' not in front_matter:
            front_matter['data'] = {}
        
        if not isinstance(front_matter['data'], dict):
            return False
        
        # Add minified code as custom fields (we'll display these in a custom template)
        if minified_code:
            front_matter['data']['minifiedCode'] = minified_code
            updated = True
        
        if minified_bundled_code:
            front_matter['data']['minifiedBundledCode'] = minified_bundled_code
            updated = True
        
        if not updated:
            return False
        
        # Re-serialize YAML front matter with literal blocks for multiline strings
        new_front_matter_str = yaml.dump(
            front_matter,
            Dumper=_LiteralDumper,
            default_flow_style=False,
            allow_unicode=True,
            sort_keys=False,
        )
        
        # Write updated content
        # Emit well-formed front matter blocks; Jekyll requires newlines after delimiters
        new_content = f'---\n{new_front_matter_str}---\n{body}'
        with open(markdown_file, 'w', encoding='utf-8') as f:
            f.write(new_content)
        
        return True
    except Exception as e:
        print(f"Error processing {markdown_file}: {e}", file=sys.stderr)
        return False


def main():
    markdown_dir = Path('_jekyll')
    minified_dir = Path('cp-algo/min')
    minified_bundled_dir = Path('.competitive-verifier/minified-bundled')
    
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
    print(f"  Markdown dir: {markdown_dir}")
    print(f"  Minified source dir: {minified_dir}")
    print(f"  Minified bundled dir: {minified_bundled_dir}")
    
    count = 0
    # Find all markdown files
    md_files = list(markdown_dir.rglob('*.md'))
    print(f"  Found {len(md_files)} markdown files")
    
    for md_file in md_files:
        # Get relative path without .md extension, then drop original source extension (.hpp/.cpp/.h)
        rel_path = md_file.relative_to(markdown_dir)
        path_without_ext = str(rel_path)[:-3]  # Remove trailing .md

        # Determine if this is a test file or library file
        is_test_file = path_without_ext.startswith('verify/')
        
        # Drop the original source extension so we don't end up with double extensions
        original_ext = None
        for src_ext in ('.hpp', '.cpp', '.h'):
            if path_without_ext.endswith(src_ext):
                path_without_ext = path_without_ext[: -len(src_ext)]
                original_ext = src_ext
                break

        # Strip cp-algo/ prefix if present since min/min-bundled dirs don't duplicate it
        # e.g., _jekyll/cp-algo/math/fft.md -> look for cp-algo/min/math/fft.hpp
        path_in_min = path_without_ext
        if path_in_min.startswith('cp-algo/'):
            path_in_min = path_in_min[8:]  # Remove 'cp-algo/' prefix
        
        minified_code = None
        minified_bundled_code = None
        
        if is_test_file:
            # For test files, generate minified versions on-the-fly from the source file
            # Reconstruct source file path
            if original_ext:
                source_file = Path(path_without_ext + original_ext)
                if source_file.exists():
                    # Read the source and minify directly
                    with open(source_file, 'r', encoding='utf-8') as f:
                        source_code = f.read()
                    minified_code = minify_code(source_code)
                    
                    # Generate bundled+minified version
                    minified_bundled_code = bundle_and_minify(source_file)
        else:
            # For library files, use pre-generated minified versions
            # Try to find corresponding minified source file
            possible_extensions = ['hpp', 'cpp', 'h']
            
            for ext in possible_extensions:
                if minified_code is None and minified_dir.exists():
                    minified_file = minified_dir / f"{path_in_min}.{ext}"
                    if minified_file.exists():
                        with open(minified_file, 'r', encoding='utf-8') as f:
                            minified_code = f.read()
                        break
            
            # Try to find corresponding minified bundled file
            # Try both with and without cp-algo/ prefix for backwards compatibility
            for ext in possible_extensions:
                if minified_bundled_code is None and minified_bundled_dir.exists():
                    # Try without prefix first (correct structure after fix)
                    minified_bundled_file = minified_bundled_dir / f"{path_in_min}.{ext}"
                    if minified_bundled_file.exists():
                        with open(minified_bundled_file, 'r', encoding='utf-8') as f:
                            minified_bundled_code = f.read()
                        break
                    # Also try with cp-algo/ prefix (old nested structure)
                    minified_bundled_file = minified_bundled_dir / f"cp-algo/{path_in_min}.{ext}"
                    if minified_bundled_file.exists():
                        with open(minified_bundled_file, 'r', encoding='utf-8') as f:
                            minified_bundled_code = f.read()
                        break
        
        # Only inject if we found at least one minified version
        if (minified_code or minified_bundled_code):
            if inject_minified_to_markdown(md_file, minified_code, minified_bundled_code):
                count += 1
                print(f"  Updated: {path_without_ext}")
            else:
                print(f"  Failed to inject into: {path_without_ext}")
        else:
            # Debug: show why no files were found
            if count == 0 or len(md_file.relative_to(markdown_dir).parts) == 3:  # Print first few
                print(f"  No minified files found for: {path_without_ext}")
                print(f"    Looked for: {path_in_min}.*")
    
    print(f"\nUpdated {count} documentation files")
    return 0


if __name__ == '__main__':
    sys.exit(main())
