# Due to the directory structure of QFlex
# go up in hierarchy twice
import sys
sys.path.insert(1, '../../')

from qflexcirq import QFlexGrid


def test_create_rectangular():
    sgrid = QFlexGrid.create_rectangular(2, 3)
    assert (sgrid == "111\n111")


def test_get_qflex_file_contents():
    sgrid = "  101   \n   110  "

    qgrid = QFlexGrid(qflex_grid_strings=sgrid)
    with open(qgrid.temp_file_if.fullpath, "r") as file:
        fcont = "".join(file.readlines()).strip()

    assert (fcont == "".join(["1 0 1\n", "1 1 0\n"]).strip())


def test_get_qubits_off():
    qub_off = QFlexGrid.get_qubits_off("010\n   101")
    assert (qub_off == [(0, 0), (0, 2), (1, 1)])
    assert (len(qub_off) == 3)
