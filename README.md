░█▀▀▀█ ░█▄─░█ ▀▀█▀▀ █▀▄▀█ ░█▀▀▀█ ░█▄─░█ ▀▀█▀▀<br>
░█──░█ ░█░█░█ ─░█── █─▀─█ ░█──░█ ░█░█░█ ─░█──<br>
░█▄▄▄█ ░█──▀█ ─░█── ▀───▀ ░█▄▄▄█ ░█──▀█ ─░█──<br>
<br>
Package for analyzing split reads from ONT (Oxford Nanopore Technologies) alignment data. 

## Install
This simple script will get things going:
```bash
pip install ontmont
```

## Test
Please check out whether the unit tests work before you begin using the module, by:
```bash
cd /path/to/ontmont
pytest
```
Ideally you should get some output like:
```bash
============================ test session starts =============================
platform linux -- Python 3.11.8, pytest-8.1.1, pluggy-1.4.0
rootdir: /data1/shahs3/users/chois7/projects/ontmont
configfile: pyproject.toml
plugins: cov-5.0.0, anyio-4.3.0
collected 49 items

tests/test_collect.py ........                                         [ 16%]
tests/test_datatypes.py ...........                                    [ 38%]
tests/test_irs.py .......................                              [ 85%]
tests/test_process.py ....                                             [ 93%]
tests/test_utils.py ...                                                [100%]

============================= 49 passed in 2.49s =============================
```

## Usage
Documentation about the modules are in the [Read the Docs](https://ontmont.readthedocs.io/en/latest/) pages.

## Development
Do build and upload to pypi:
    1. `python -m build`
    2. `python -m twine upload --repository pypi dist/*`
