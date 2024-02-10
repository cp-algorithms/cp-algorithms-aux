---
data:
  _extendedDependsOn: []
  _extendedRequiredBy: []
  _extendedVerifiedWith: []
  _isVerificationFailed: false
  _pathExtension: py
  _verificationStatusIcon: ':warning:'
  attributes:
    links: []
  bundledCode: "Traceback (most recent call last):\n  File \"/home/runner/.local/lib/python3.10/site-packages/onlinejudge_verify/documentation/build.py\"\
    , line 71, in _render_source_code_stat\n    bundled_code = language.bundle(stat.path,\
    \ basedir=basedir, options={'include_paths': [basedir]}).decode()\n  File \"/home/runner/.local/lib/python3.10/site-packages/onlinejudge_verify/languages/python.py\"\
    , line 96, in bundle\n    raise NotImplementedError\nNotImplementedError\n"
  code: "import networkx as nx\nimport matplotlib.pyplot as plt\nfor seed in range(100):\n\
    \    G = nx.gnp_random_graph(6, 0.6, seed=seed)\n\n    plt.subplot(111)\n    nx.draw(G)\
    \   # default spring_layout\n    print(seed)\n    plt.show()\n"
  dependsOn: []
  isVerificationFile: false
  path: images/second_best_mst/random_graph.py
  requiredBy: []
  timestamp: '1970-01-01 00:00:00+00:00'
  verificationStatus: LIBRARY_NO_TESTS
  verifiedWith: []
documentation_of: images/second_best_mst/random_graph.py
layout: document
redirect_from:
- /library/images/second_best_mst/random_graph.py
- /library/images/second_best_mst/random_graph.py.html
title: images/second_best_mst/random_graph.py
---
