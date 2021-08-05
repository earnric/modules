#! /bin/env python3

import sys

from kepler.code.main import Kepler, KeplerData
from kepler.cmdloop import KeplerCmdLoop

class KeplerShell(KeplerCmdLoop):
    def __init__(self, *args):
        k = Kepler(*args)
        d = KeplerData()
        super().__init__(k, d)

def run(*args):
    KeplerShell(*args).cmdloop()

if __name__ == "__main__":
    run(*sys.argv[1:])
