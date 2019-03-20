import numpy as np
import itertools as it

filename = '../circuits/ben_11_24_0.txt'

#ordered_sites = [  0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11,
#                  12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23,
#                  24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35,
#                  36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47,
#                  48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59,
#                  60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71,
#                  72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83,
#                  84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95,
#                  96, 97, 98, 99,100,101,102,103,104,105,106,107,
#                 108,109,110,111,112,113,114,115,116,117,118,119,
#                 120,121,122,123,124,125,126,127,128,129,130,131]
#ordered_sites = [                      5,  6,                    
#                                  16, 17, 18, 19,                
#                                  28, 29, 30, 31, 32,            
#                          38, 39, 40, 41, 42, 43, 44, 45,        
#                      49, 50, 51, 52, 53, 54, 55, 56, 57,        
#                      61, 62, 63, 64, 65, 66, 67, 68,            
#                      73, 74, 75, 76, 77, 78, 79,                
#                          86, 87, 88, 89, 90,                    
#                              99,100,101,                        
#                                                                 
#                                                                ]
#ordered_sites = [                                                
#                                                                 
#                                                                 
#                                                                 
#                                                                 
#                                                                 
#                                                                 
#                                  88, 89, 90                     
#                                                                 
#                                                                 
#                                                                ]
ordered_sites = [
                
               
              
             
                                      
                                  76, 77, 78,     
                                      89    
                
                 
                                                                ]
sites_map = {s: i for i, s in enumerate(ordered_sites)}

# Ignore all gates involving qubits not in ordered_sites. Ordered sites will
# map onto {0..63}

with open(filename, 'r') as f:
  lines = f.readlines()

print(len(ordered_sites))

for line in lines[1:]:
  l = line.strip().split(' ')
  str_print = ''
  if len(l)==3:
    if int(l[2]) in ordered_sites:
      str_print += l[0]
      str_print += ' '
      str_print += l[1]
      str_print += ' '
      str_print += str(sites_map[int(l[2])])
      print(str_print)
  if len(l)==4:
    if (int(l[2]) in ordered_sites) and (int(l[3]) in ordered_sites):
      str_print += l[0]
      str_print += ' '
      str_print += l[1]
      str_print += ' '
      str_print += str(sites_map[int(l[2])])
      str_print += ' '
      str_print += str(sites_map[int(l[3])])
      print(str_print)

