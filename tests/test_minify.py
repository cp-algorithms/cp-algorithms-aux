import importlib.util
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest


spec = importlib.util.spec_from_file_location(
    'minify', Path(__file__).resolve().parents[1] / 'cp-algo/minify.py')
minifier = importlib.util.module_from_spec(spec)
spec.loader.exec_module(minifier)


class MinifierTests(unittest.TestCase):
    def test_existing_compound_operators(self):
        code = 'x >>= 1; y += 2; ++z; p->member; a <=> b;'
        result = minifier.minify_cpp(code)
        for token in ['>>=', '+=', '++', '->', '<=>']:
            self.assertIn(token, result)

    def test_string_whitespace(self):
        literal = '"a + + b > = c // text"'
        self.assertIn(literal, minifier.minify_cpp('auto text = ' + literal + ';'))

    def test_compiled_token_boundaries(self):
        compiler = shutil.which(os.environ.get('CXX', 'g++'))
        if not compiler:
            self.skipTest('C++ compiler unavailable')
        # The specialization reproduces the minified Strassen-header failure.
        source = r'''
#include <cassert>
template<class> constexpr bool marked = false;
template<int n> struct value {};
template<class T> struct row {};
template<int n> constexpr bool marked<row<value<n>>> = n > 0;
int main() {
    static_assert(marked<row<value<1>>>);
    int a = 3, b = 2;
    assert(a + + b == 5 && a - - b == 5);
    assert(a == 3 && b == 2);
    int x = a +
        + b;
    unsigned
        long y = 8;
    int/**/z = 4;
    y >>= 1;
    assert(x == 5 && y == 4 && z == 4);
}
'''
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for name, text in [('original', source),
                               ('minified', minifier.minify_cpp(source))]:
                path = root / (name + '.cpp')
                path.write_text(text)
                executable = root / name
                subprocess.run([compiler, '-std=c++23', str(path), '-o', str(executable)],
                               check=True, capture_output=True, text=True)
                subprocess.run([str(executable)], check=True)


if __name__ == '__main__':
    unittest.main()
