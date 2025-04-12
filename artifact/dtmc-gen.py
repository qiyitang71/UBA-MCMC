import random
from decimal import Decimal
import sys

def generate_prism_lmc(num_states=1000, num_labels=4):
    labels = ["sigma", "pi", "hash", "dollar"]
    prism_model = "dtmc\n\nmodule random\n"
    prism_model += "\t// local state\n"
    prism_model += f"\ts : [0..{num_states - 1}] init 0;\n\n"

    for s in range(num_states):
        # Group states by their labels
        states_by_label = {label: [] for label in range(num_labels)}
        for state in range(num_states):
            label = state % num_labels
            states_by_label[label].append(state)

        # Randomly select a random number of targets (at most 4)
        num_transitions = random.randint(1, min(num_labels, 4))  # Up to 4 transitions
        available_labels = [label for label in range(num_labels)]

        # Randomly select unique target states from different labels
        target_states = []
        while len(target_states) < num_transitions:
            label = random.choice(available_labels)
            if states_by_label[label]:  # Check if there are available states with this label
                target = random.choice(states_by_label[label])
                if target not in target_states:
                    target_states.append(target)
                    available_labels.remove(label)  # Remove used label

        # Normalize probabilities
        probabilities = [Decimal(1) / len(target_states) for _ in range(len(target_states))]
        total = sum(probabilities)
        probabilities = [p / total for p in probabilities]

        # Ensure probabilities sum to exactly 1
        probabilities[-1] += 1 - sum(probabilities)

        # Construct the transition line
        transitions = " + ".join(
            f"{prob:.5f} : (s'={target})"
            for prob, target in zip(probabilities, target_states)
        )
        prism_model += f"\t[] s={s} -> {transitions};\n"

    prism_model += "endmodule\n\n"

    # Define labels
    prism_model += "// Labels\n"
    for i, label in enumerate(labels):
        prism_model += f"label \"{label}\" = mod(s,{num_labels})={i};\n"

    return prism_model

if __name__ == "__main__":
    # Check if a file name is provided as argument
    if len(sys.argv) != 2:
        print("Usage: python3 generate_lmc.py <output_file>")
        sys.exit(1)

    output_file = sys.argv[1]
    
    # Generate a random labelled Markov chain with default parameters
    model = generate_prism_lmc()
    
    # Save the model to the specified file
    with open(output_file, "w") as f:
        f.write(model)
    
    print(f"Random Markov chain model with distinct labels for transitions generated and saved to '{output_file}'.")