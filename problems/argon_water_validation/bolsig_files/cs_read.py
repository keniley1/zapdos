import numpy as np
import matplotlib.pyplot as plt

# Ars -> Arss
filename = 'Ars_to_Arss.txt'

in_reaction = False
reading_data = False
with open(filename, 'r') as cs_file:
    for line in cs_file:
        if ("EXCITATION" in line):
            in_reaction = True
            continue

        if ("----" in line and in_reaction and not reading_data):
            reading_data = True
            continue
        elif ("----" in line and in_reaction and reading_data):
            in_reaction = False
            reading_data = False

        if (reading_data):
            print(line)


            
