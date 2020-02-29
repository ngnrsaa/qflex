import sys
sys.path.insert(1, '../../')

from qflexcirq import qflex_cirq_example as qce


def test_example():
    # change CWD - the hacky way
    import os
    init_cwd = os.getcwd()
    os.chdir("../..")

    # If it does not crash, then it is fine
    qce.main()

    # Return CWD to its starting state
    os.chdir(init_cwd)
