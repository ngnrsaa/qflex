import qflex

with open("./circuits/ben_11_16_0.txt", "r") as f:
    circuit_content = f.readlines()

with open("./ordering/bristlecone_48.txt", "r") as f:
    ordering_content = f.readlines()

with open("./grid/bristlecone_48.txt", "r") as f:
    grid_content = f.readlines()

amplitudes = qflex.simulate(circuit_content, ordering_content, grid_content)

# hard coded
input_initial_state = "XXXX"

print("This is from Python")
for amp in amplitudes:
    state = amp[0]
    amplitude = complex(amp[1])

    print(input_initial_state + " --> " + state + ": " + \
            str(amplitude.real) + " " + str(amplitude.imag))