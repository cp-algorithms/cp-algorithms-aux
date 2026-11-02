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

Two workflows use these scripts:

1. `.github/workflows/verify.yml` bundles the files with `oj-resolve`, runs `generate_minified.py` and runs `inject_minified_docs.py`.
2. After a successful `verify` run on `main`, `.github/workflows/publish.yml` runs the scripts again. If `cp-algo/min/` changed, it commits the update and skips the site build. Otherwise it builds and deploys the Jekyll site.

The bot's commit does not start a new workflow run. After a minified update, the site is rebuilt by the next push to `main`.

## Website Display

Each file on the [library site](https://lib.cp-algorithms.com) shows up to 4 code tabs:
- **default**: original source code
- **bundled**: all dependencies inlined
- **minified**: minified source, about 35% smaller in total than the original
- **minified + bundled**: the bundled version, minified

## Local Usage

After cloning the repository, you can directly use minified versions from the `cp-algo/min/` directory:

```bash
# Copy minified version of a file
cp cp-algo/min/graph/mst.hpp your_solution.hpp

# Or use it inline in competitive programming
#include "cp-algo/min/graph/mst.hpp"
```
