import random
import sys
from decimal import Decimal, getcontext, ROUND_DOWN

# Set precision high enough
getcontext().prec = 28

def generate_prism_lmc(num_states=20, num_labels=4):
    labels = ["sigma", "pi", "hash", "dollar"]
    prism_model = "dtmc\n\nmodule random\n"
    prism_model += "\t// local state\n"
    prism_model += f"\ts : [0..{num_states - 1}] init 0;\n\n"

    for s in range(num_states):
        # Group states by label residue mod num_labels
        states_by_label = {label: [] for label in range(num_labels)}
        for state in range(num_states):
            label = state % num_labels
            states_by_label[label].append(state)

        # Pick one state from each label residue
        target_states = []
        used_states = set()
        for label in range(num_labels):
            candidates = states_by_label[label]
            if not candidates:
                raise ValueError(f"No state with label mod {num_labels} == {label}")
            chosen = random.choice(candidates)
            target_states.append(chosen)
            used_states.add(chosen)

        # Add additional unique states (optional)
        extra_states = [
            state for state in range(num_states) if state not in used_states
        ]
        extra_targets = random.sample(extra_states, random.randint(0, len(extra_states)))
        target_states.extend(extra_targets)

        n = len(target_states)
        base_prob = Decimal('1') / Decimal(n)

        # Round down all but last
        raw_probs = [base_prob] * n
        fixed_strs = [p.quantize(Decimal('1.00000'), rounding=ROUND_DOWN) for p in raw_probs[:-1]]
        sum_fixed = sum(fixed_strs)
        final_prob = Decimal('1.00000') - sum_fixed
        fixed_strs.append(final_prob.quantize(Decimal('1.00000')))

        # Convert to strings for printing
        transition_strs = [
            f"{p} : (s'={t})" for p, t in zip(fixed_strs, target_states)
        ]
        transitions = " + ".join(transition_strs)
        prism_model += f"\t[] s={s} -> {transitions};\n"

    prism_model += "endmodule\n\n"

    # Add labels
    for i, label in enumerate(labels):
        prism_model += f"label \"{label}\" = mod(s,{num_labels})={i};\n"

    return prism_model

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python3 generate_lmc.py <output_file>")
        sys.exit(1)

    output_file = sys.argv[1]
    model = generate_prism_lmc()

    with open(output_file, "w") as f:
        f.write(model)
