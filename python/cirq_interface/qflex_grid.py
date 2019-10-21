import tempfile
import os

class QFlexGrid():
    BRISTLECONE48 = """000001100000
                    000011110000
                    000111111000
                    001111111100
                    001111111100
                    001111111100
                    000111111000
                    000011110000
                    000001100000
                    000000000000
                    000000000000"""#11 lines of 12 cols

    BRISTLECONE60 = """000001100000
                    000011110000
                    000111111000
                    001111111100
                    011111111110
                    011111111110
                    001111111100
                    000111111000
                    000011110000
                    000001100000
                    000000000000"""

    BRISTLECONE70 = """000001100000
                    000011110000
                    000111111000
                    001111111100
                    011111111110
                    011111111110
                    011111111110
                    001111111100
                    000111111000
                    000011110000
                    000001100000"""#11 lines of 12 cols

    def __init__(self, qflex_grid_strings = BRISTLECONE70):
        # Behind the scene, this class creates a temporary file for each object
        self._file_handle = tempfile.mkstemp()

        with open(self._file_handle[1], "w") as f:
            # I do have the file handle anyway...
            print(qflex_grid_strings, file = f)

    def __del__(self):
        # The destructor removes the temporary file

        # if open, close the file handle
        try:
            os.close(self._file_handle[0])
        except OSError as e:
            if e.errno == 9:
                # if it was closed before
                pass
            else:
                raise e


        # remove the temporary file from disk
        os.remove(self._file_handle[1])

    @staticmethod
    def create_rectangular(sizex, sizey):
        regular = ""

        for x in range(sizex):
            line = "1" * sizey

            if x > 0 :
                regular += "\n"

            regular += line

        return regular

    # @staticmethod
    # def get_qflex_file_contents(qflex_grid_string):
    #     gdata = qflex_grid_string.replace("0", "0 ") \
    #         .replace("1", "1 ")
    #
    #     grid_data = [x.strip() + "\n" for x in gdata.split("\n")]
    #
    #     return grid_data

    @staticmethod
    def get_qubits_off(qflex_grid_string):
        qubits_off = []

        for i, x in enumerate(qflex_grid_string.split("\n")):
            for j, y in enumerate(x.strip()):
                if y == "0":
                    qubits_off.append((i, j))

        return qubits_off