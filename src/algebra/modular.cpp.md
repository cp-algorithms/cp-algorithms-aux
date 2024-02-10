---
data:
  _extendedDependsOn: []
  _extendedRequiredBy: []
  _extendedVerifiedWith: []
  _isVerificationFailed: false
  _pathExtension: cpp
  _verificationStatusIcon: ':warning:'
  attributes:
    links:
    - https://en.wikipedia.org/wiki/Berlekamp-Rabin_algorithm
  bundledCode: "#line 1 \"src/algebra/modular.cpp\"\nnamespace algebra { // modular\n\
    \    template<int m>\n    struct modular {\n        // https://en.wikipedia.org/wiki/Berlekamp-Rabin_algorithm\n\
    \        // solves x^2 = y (mod m) assuming m is prime in O(log m).\n        //\
    \ returns nullopt if no sol.\n        optional<modular> sqrt() const {\n     \
    \       static modular y;\n            y = *this;\n            if(r == 0) {\n\
    \                return 0;\n            } else if(bpow(y, (m - 1) / 2) != modular(1))\
    \ {\n                return nullopt;\n            } else {\n                while(true)\
    \ {\n                    modular z = rng();\n                    if(z * z == *this)\
    \ {\n                        return z;\n                    }\n              \
    \      struct lin {\n                        modular a, b;\n                 \
    \       lin(modular a, modular b): a(a), b(b) {}\n                        lin(modular\
    \ a): a(a), b(0) {}\n                        lin operator * (const lin& t) {\n\
    \                            return {\n                                a * t.a\
    \ + b * t.b * y,\n                                a * t.b + b * t.a\n        \
    \                    };\n                        }\n                    } x(z,\
    \ 1); // z + x\n                    x = bpow(x, (m - 1) / 2);\n              \
    \      if(x.b != modular(0)) {\n                        return x.b.inv();\n  \
    \                  }\n                }\n            }\n        }\n        \n\
    \        uint64_t r;\n        constexpr modular(): r(0) {}\n        constexpr\
    \ modular(int64_t rr): r(rr % m) {r = min<uint64_t>(r, r + m);}\n        modular\
    \ inv() const {return bpow(*this, m - 2);}\n        modular operator - () const\
    \ {return min(-r, m - r);}\n        modular operator * (const modular &t) const\
    \ {return r * t.r;}\n        modular operator / (const modular &t) const {return\
    \ *this * t.inv();}\n        modular& operator += (const modular &t) {r += t.r;\
    \ r = min<uint64_t>(r, r - m); return *this;}\n        modular& operator -= (const\
    \ modular &t) {r -= t.r; r = min<uint64_t>(r, r + m); return *this;}\n       \
    \ modular operator + (const modular &t) const {return modular(*this) += t;}\n\
    \        modular operator - (const modular &t) const {return modular(*this) -=\
    \ t;}\n        modular& operator *= (const modular &t) {return *this = *this *\
    \ t;}\n        modular& operator /= (const modular &t) {return *this = *this /\
    \ t;}\n        \n        auto operator <=> (const modular &t) const = default;\n\
    \        \n        explicit operator int() const {return r;}\n        int64_t\
    \ rem() const {return 2 * r > m ? r - m : r;}\n\n        static constexpr uint64_t\
    \ mm = (uint64_t)m * m;\n        void add_unsafe(uint64_t t) {r += t; r = min<uint64_t>(r,\
    \ r - mm);}\n        modular& normalize() {if(r >= m) r %= m; return *this;}\n\
    \    };\n    \n    template<int m>\n    istream& operator >> (istream &in, modular<m>\
    \ &x) {\n        return in >> x.r;\n    }\n    \n    template<int m>\n    ostream&\
    \ operator << (ostream &out, modular<m> const& x) {\n        return out << x.r\
    \ % m;\n    }\n}\n"
  code: "namespace algebra { // modular\n    template<int m>\n    struct modular {\n\
    \        // https://en.wikipedia.org/wiki/Berlekamp-Rabin_algorithm\n        //\
    \ solves x^2 = y (mod m) assuming m is prime in O(log m).\n        // returns\
    \ nullopt if no sol.\n        optional<modular> sqrt() const {\n            static\
    \ modular y;\n            y = *this;\n            if(r == 0) {\n             \
    \   return 0;\n            } else if(bpow(y, (m - 1) / 2) != modular(1)) {\n \
    \               return nullopt;\n            } else {\n                while(true)\
    \ {\n                    modular z = rng();\n                    if(z * z == *this)\
    \ {\n                        return z;\n                    }\n              \
    \      struct lin {\n                        modular a, b;\n                 \
    \       lin(modular a, modular b): a(a), b(b) {}\n                        lin(modular\
    \ a): a(a), b(0) {}\n                        lin operator * (const lin& t) {\n\
    \                            return {\n                                a * t.a\
    \ + b * t.b * y,\n                                a * t.b + b * t.a\n        \
    \                    };\n                        }\n                    } x(z,\
    \ 1); // z + x\n                    x = bpow(x, (m - 1) / 2);\n              \
    \      if(x.b != modular(0)) {\n                        return x.b.inv();\n  \
    \                  }\n                }\n            }\n        }\n        \n\
    \        uint64_t r;\n        constexpr modular(): r(0) {}\n        constexpr\
    \ modular(int64_t rr): r(rr % m) {r = min<uint64_t>(r, r + m);}\n        modular\
    \ inv() const {return bpow(*this, m - 2);}\n        modular operator - () const\
    \ {return min(-r, m - r);}\n        modular operator * (const modular &t) const\
    \ {return r * t.r;}\n        modular operator / (const modular &t) const {return\
    \ *this * t.inv();}\n        modular& operator += (const modular &t) {r += t.r;\
    \ r = min<uint64_t>(r, r - m); return *this;}\n        modular& operator -= (const\
    \ modular &t) {r -= t.r; r = min<uint64_t>(r, r + m); return *this;}\n       \
    \ modular operator + (const modular &t) const {return modular(*this) += t;}\n\
    \        modular operator - (const modular &t) const {return modular(*this) -=\
    \ t;}\n        modular& operator *= (const modular &t) {return *this = *this *\
    \ t;}\n        modular& operator /= (const modular &t) {return *this = *this /\
    \ t;}\n        \n        auto operator <=> (const modular &t) const = default;\n\
    \        \n        explicit operator int() const {return r;}\n        int64_t\
    \ rem() const {return 2 * r > m ? r - m : r;}\n\n        static constexpr uint64_t\
    \ mm = (uint64_t)m * m;\n        void add_unsafe(uint64_t t) {r += t; r = min<uint64_t>(r,\
    \ r - mm);}\n        modular& normalize() {if(r >= m) r %= m; return *this;}\n\
    \    };\n    \n    template<int m>\n    istream& operator >> (istream &in, modular<m>\
    \ &x) {\n        return in >> x.r;\n    }\n    \n    template<int m>\n    ostream&\
    \ operator << (ostream &out, modular<m> const& x) {\n        return out << x.r\
    \ % m;\n    }\n}\n"
  dependsOn: []
  isVerificationFile: false
  path: src/algebra/modular.cpp
  requiredBy: []
  timestamp: '2024-02-10 16:40:11+01:00'
  verificationStatus: LIBRARY_NO_TESTS
  verifiedWith: []
documentation_of: src/algebra/modular.cpp
layout: document
redirect_from:
- /library/src/algebra/modular.cpp
- /library/src/algebra/modular.cpp.html
title: src/algebra/modular.cpp
---
