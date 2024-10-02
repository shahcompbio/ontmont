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
======================= test session starts =======================
platform linux -- Python 3.11.8, pytest-8.1.1, pluggy-1.4.0
rootdir: /data1/shahs3/users/chois7/projects/ontmont
configfile: pyproject.toml
plugins: cov-5.0.0, anyio-4.3.0
collected 46 items

tests/test_collect.py .....                                 [ 10%]
tests/test_datatypes.py ...........                         [ 34%]
tests/test_irs.py .......................                   [ 84%]
tests/test_process.py ....                                  [ 93%]
tests/test_utils.py ...                                     [100%]

======================= 46 passed in 2.50s ========================
```

## Usage
Documentation about the modules are in the [Read the Docs](https://ontmont.readthedocs.io/en/latest/) pages.


