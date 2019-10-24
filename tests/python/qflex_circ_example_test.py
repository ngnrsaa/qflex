import sys
sys.path.insert(1, '../../')

from python import qflex_cirq_example as qce

def test_example():
    # change CWD - the hacky way
    import os
    os.chdir("../..")

    # If it does not crash, then it is fine
    qce.main()