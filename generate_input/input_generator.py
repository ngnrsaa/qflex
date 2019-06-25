import numpy as np

I = 10
J = 10
K = 2
fidelity = 0.005
filename = 'circuits/circuit_q53c0d4.qflex'
num_qubits = 53
num_qubits_A = 6
num_qubits_B = num_qubits - num_qubits_A
num_strings = 1000
#num_strings_A = 2**num_qubits_A
random_seed = 1

################################   DONE!   ####################################

np.random.seed(random_seed)

input_string = ''
for i in range(num_qubits):
  input_string += '0'

common_string = '{} {} {} {} {}'.format(I, J, K, fidelity, filename)

strings_A = []
for i in range(2**num_qubits_A):
  s = ''
  for n in range(num_qubits_A):
    c = str(int(i&(1<<n)>0))
    s = s + c
  strings_A.append(s)

strings_B = {}
strings_left = num_strings
while strings_left:
  s = ''
  for n in range(num_qubits_B):
    s = s + str(int(np.random.randint(2)))
  if s not in strings_B:
    strings_B[s] = 0 # dummy value
    strings_left -= 1

strings_B_list = [s for s in strings_B]

print(2 + len(strings_A))
print(common_string)
for i in range(num_strings):
  out_string = str(input_string) + ' '
  out_string += strings_B_list[i] + ' '
  for j in range(len(strings_A)):
    out_string += strings_A[j] + ' '
  print(out_string)
