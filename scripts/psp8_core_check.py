"""The check of psp8_core_test.sh: the psp8 file of the output NEW.out has
the model core (fchrg > 0 and the model core block, of the shape of that of
BASE.out, the unchanged input, whose model core is written).

usage: psp8_core_check.py BASE.out NEW.out
"""

import sys

from pyoncvpsp.io import OncvpspTextParser, parse_psp8

base, new = (
    parse_psp8(OncvpspTextParser(path=path, check=False).psp8_string)
    for path in sys.argv[1:]
)
if base["model_core"] is None:
    sys.exit(
        "psp8_core_check.py: the psp8 file of the unchanged input has no model core"
    )
if new["fchrg"] <= 0:
    sys.exit(f"psp8_core_check.py: fchrg={new['fchrg']}, no model core")
if new["model_core"] is None or new["model_core"].shape != base["model_core"].shape:
    sys.exit("psp8_core_check.py: the psp8 file has no model core block")
print(f"the psp8 file has the model core, fchrg={new['fchrg']}")
