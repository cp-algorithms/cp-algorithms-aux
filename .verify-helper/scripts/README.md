# Minified Versions Setup

This directory contains scripts to generate and manage minified code versions for the library.

## Files

- `generate_minified.py`: Generates minified versions from source files
- `inject_minified_docs.py`: Injects minified code into documentation markdown files

## How it works

1. **generate_minified.py**: Reads all files from `cp-algo/` and creates minified versions in:
   - `.competitive-verifier/minified/cp-algo/` (for CI/documentation)
   - `cp-algo/min/` (committed to repo for direct access)
2. **inject_minified_docs.py**: Takes the minified files from `cp-algo/min/` and adds them to the documentation front matter as `minifiedCode` fields
3. The Jekyll template in `_includes/document_body.html` displays the minified code when available

## GitHub Actions Integration

The workflow in `.github/workflows/verify.yml` automatically:
1. Generates bundled files (via `oj-resolve`)
2. Calls `generate_minified.py` to create minified versions in both locations
3. Calls `inject_minified_docs.py` to add them to documentation
4. Commits and pushes the `cp-algo/min/` directory updates to the repository
5. Builds the Jekyll site with all three code versions displayed

## Website Display

On the library website (https://lib.cp-algorithms.com), each file will now show three code tabs:
- **default**: Original source code
- **bundled**: Bundled version with all dependencies inlined
- **minified**: Minified version (typically 35-50% size reduction from source)

## Local Usage

After cloning the repository, you can directly use minified versions from the `cp-algo/min/` directory:

```bash
# Copy minified version of a file
cp cp-algo/min/graph/mst.hpp your_solution.hpp

# Or use it inline in competitive programming
#include "cp-algo/min/graph/mst.hpp"
```
