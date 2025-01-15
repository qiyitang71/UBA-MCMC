#!/usr/bin/env python3

labels = ["sigma", "pi", "hash", "dollar"]

# Open the file in write mode
with open("labels.txt", "w") as file:
    for i in range(500):  # Two sets: _0 and _1
        for j, label in enumerate(labels):
            file.write(f'label "{label}_{i}" = s={i * len(labels) + j};\n')