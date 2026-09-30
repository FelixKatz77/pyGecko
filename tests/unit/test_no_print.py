'''The library never prints.

Progress goes to logging, advisories to warnings and errors are raised, so a caller (a notebook,
a server, an agent) decides what reaches a console. Checked on the syntax tree so a print call in
any module fails here, whether or not a test happens to execute it.
'''

import ast
from pathlib import Path

import pytest

import pygecko

PACKAGE = Path(pygecko.__file__).parent
SOURCES = sorted(PACKAGE.rglob('*.py'))


@pytest.mark.parametrize('path', SOURCES, ids=lambda path: str(path.relative_to(PACKAGE)))
def test_module_contains_no_print_call(path):
    tree = ast.parse(path.read_text(encoding='utf-8'))
    lines = [node.lineno for node in ast.walk(tree)
             if isinstance(node, ast.Call) and isinstance(node.func, ast.Name) and node.func.id == 'print']
    assert lines == [], f'print() at lines {lines}'
