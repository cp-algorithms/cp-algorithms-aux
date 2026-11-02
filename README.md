# e-maxx-eng-aux

[![Actions Status](https://github.com/cp-algorithms/cp-algorithms-aux/workflows/verify/badge.svg)](https://github.com/cp-algorithms/cp-algorithms-aux/actions)
[![GitHub Pages](https://img.shields.io/static/v1?label=GitHub+Pages&message=CP-Algorithms+library+&color=brightgreen&logo=github)](https://lib.cp-algorithms.com/)

Auxiliary stuff for article sources.

## Documentation previews

Successful verification builds a `pr-preview` artifact for pull requests. Only
open PRs from branches in this repository are published at
`https://lib.cp-algorithms.com/<PR number>/`; the publisher checks the current
head commit before publishing. Fork previews remain downloadable artifacts,
even after a maintainer approves their CI run. Untrusted fork HTML must not be
served from the production site's origin. Closing a PR removes its published
preview. Publishing fork previews would require a separate, isolated origin.
