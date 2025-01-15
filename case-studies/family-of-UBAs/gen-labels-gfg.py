#!/usr/bin/env python3

labels = ["sigma", "pi", "hash", "dollar"]

result = []
for i in range(500):  # Two sets: _0 and _1
    for label in labels:
        result.append(f'"{label}_0" <- "{label}_{i}"')

# Join the result with commas and print
# Write the result to a file
with open("mappings-gfg.txt", "w") as file:
    file.write(", ".join(result))
