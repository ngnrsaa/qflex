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

    @staticmethod
    def create_rectangular(sizex, sizey):
        regular = ""

        for x in range(sizex):
            line = "1" * sizey

            if x >0:
                regular += "\n"

            regular += line

        return regular


    @staticmethod
    def get_qflex_file_contents(grid_string):
        gdata = grid_string.replace("0", "0 ") \
            .replace("1", "1 ")

        grid_data = [x.strip() + "\n" for x in gdata.split("\n")]

        return grid_data

    @staticmethod
    def get_qubits_off(grid_string):
        qubits_off = []

        for i, x in enumerate(grid_string.split("\n")):
            for j, y in enumerate(x.strip()):
                if y == "0":
                    qubits_off.append((i, j))

        return qubits_off