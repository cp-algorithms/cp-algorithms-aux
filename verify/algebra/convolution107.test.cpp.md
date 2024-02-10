---
data:
  _extendedDependsOn:
  - icon: ':x:'
    path: cp-algo/algebra/common.hpp
    title: cp-algo/algebra/common.hpp
  - icon: ':x:'
    path: cp-algo/algebra/fft.hpp
    title: cp-algo/algebra/fft.hpp
  - icon: ':x:'
    path: cp-algo/algebra/modular.hpp
    title: cp-algo/algebra/modular.hpp
  - icon: ':x:'
    path: cp-algo/algebra/polynomial.hpp
    title: cp-algo/algebra/polynomial.hpp
  _extendedRequiredBy: []
  _extendedVerifiedWith: []
  _isVerificationFailed: true
  _pathExtension: cpp
  _verificationStatusIcon: ':x:'
  attributes:
    '*NOT_SPECIAL_COMMENTS*': ''
    PROBLEM: https://judge.yosupo.jp/problem/convolution_mod_1000000007
    links:
    - https://judge.yosupo.jp/problem/convolution_mod_1000000007
  bundledCode: "Traceback (most recent call last):\n  File \"/home/runner/.local/lib/python3.10/site-packages/onlinejudge_verify/documentation/build.py\"\
    , line 71, in _render_source_code_stat\n    bundled_code = language.bundle(stat.path,\
    \ basedir=basedir, options={'include_paths': [basedir]}).decode()\n  File \"/home/runner/.local/lib/python3.10/site-packages/onlinejudge_verify/languages/cplusplus.py\"\
    , line 187, in bundle\n    bundler.update(path)\n  File \"/home/runner/.local/lib/python3.10/site-packages/onlinejudge_verify/languages/cplusplus_bundle.py\"\
    , line 401, in update\n    self.update(self._resolve(pathlib.Path(included), included_from=path))\n\
    \  File \"/home/runner/.local/lib/python3.10/site-packages/onlinejudge_verify/languages/cplusplus_bundle.py\"\
    , line 401, in update\n    self.update(self._resolve(pathlib.Path(included), included_from=path))\n\
    \  File \"/home/runner/.local/lib/python3.10/site-packages/onlinejudge_verify/languages/cplusplus_bundle.py\"\
    , line 401, in update\n    self.update(self._resolve(pathlib.Path(included), included_from=path))\n\
    \  File \"/home/runner/.local/lib/python3.10/site-packages/onlinejudge_verify/languages/cplusplus_bundle.py\"\
    , line 260, in _resolve\n    raise BundleErrorAt(path, -1, \"no such header\"\
    )\nonlinejudge_verify.languages.cplusplus_bundle.BundleErrorAt: cassert: line\
    \ -1: no such header\n"
  code: "#define PROBLEM \"https://judge.yosupo.jp/problem/convolution_mod_1000000007\"\
    \n#pragma GCC optimize(\"Ofast,unroll-loops\")\n#pragma GCC target(\"avx2,tune=native\"\
    )\n#include \"cp-algo/algebra/polynomial.hpp\"\n#include <bits/stdc++.h>\n\nusing\
    \ namespace std;\nusing namespace algebra;\n\nconst int mod = 1e9 + 7;\ntypedef\
    \ modular<mod> base;\ntypedef poly<base> polyn;\n\nvoid solve() {\n    int n,\
    \ m;\n    cin >> n >> m;\n    vector<base> a(n), b(m);\n    copy_n(istream_iterator<base>(cin),\
    \ n, begin(a));\n    copy_n(istream_iterator<base>(cin), m, begin(b));\n    (polyn(a)\
    \ * polyn(b)).print(n + m - 1);\n}\n\nsigned main() {\n    //freopen(\"input.txt\"\
    , \"r\", stdin);\n    ios::sync_with_stdio(0);\n    cin.tie(0);\n    int t;\n\
    \    t = 1;// cin >> t;\n    while(t--) {\n        solve();\n    }\n}\n"
  dependsOn:
  - cp-algo/algebra/polynomial.hpp
  - cp-algo/algebra/common.hpp
  - cp-algo/algebra/modular.hpp
  - cp-algo/algebra/fft.hpp
  isVerificationFile: true
  path: verify/algebra/convolution107.test.cpp
  requiredBy: []
  timestamp: '2024-02-10 18:45:35+01:00'
  verificationStatus: TEST_WRONG_ANSWER
  verifiedWith: []
documentation_of: verify/algebra/convolution107.test.cpp
layout: document
redirect_from:
- /verify/verify/algebra/convolution107.test.cpp
- /verify/verify/algebra/convolution107.test.cpp.html
title: verify/algebra/convolution107.test.cpp
---
