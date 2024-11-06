#!/usr/bin/env python3
import sys
"""
Alphabet = {0,1,2,3} 
0: sigma
1 : pi
2 : hash
3 : dollar

N is the no of boxes.
Each state of the automata is a number in range (1,2^N-2) where the binary representation of this number
is a binary string of length n, and the value of bit in i-th position  stores if q_i is in the state. 
A value of 1 denotes that q_i is in the state. For eg : for n = 3, 011
denotes the state {q_0,q_1}.  The state 0 represents the state q_$. 

Accepting transition are all labelled with 5. 
"""

def findSuccessor(state, letter, N): #gives list of successor states on reading letter from state
	if letter == 0: #sigma
		assert(state != 0)
		lastTwoBit = state & 3
		if lastTwoBit in [0,3]:
			return([0,state])
		else:
			return([0,state ^ 3])

	elif letter == 1: # pi
		assert(state != 0)
		return([0,(state & 1)*(2**(N-1))|(state >> 1)])

	elif letter ==2: # hash
		assert(state != 0)
		if state == 1:
			return([0,2**N - 2 ])
		else:
			return([0,state - (state & 1)])
	else : # dollar
		assert(state == 0)
		return([0,1])

def ubaGenerator(N):
	with open(f"uba-family-{N}.hoa",'w') as file:
		file.write(f"HOA: v1\n")
		file.write(f"name: uba-family-{N}-trans-based\n")
		NO_OF_STATES = 2**N -1
		file.write(f"States: {NO_OF_STATES}\n")
		file.write(f"Start: 1\n")
		file.write(f"acc-name: Buchi\n")
		file.write(f"Acceptance: 1 Inf(5)\n") #each acc transition is given label 5
		file.write(f"AP: 4 \"sigma\" \"pi\" \"hash\" \"dollar\"\n")
		file.write("properties: trans-labels explicit-labels trans-acc unambiguous\n")
		file.write(f"--BODY--\n")
		for state in range(NO_OF_STATES):
			file.write(f"State: {state}\n")
			if state == 0:
				file.write(f"  [3] 0 {5}\n")
				file.write(f"  [3] 1 {5}\n")
			else:
				for letter in range(3):
					for successor in findSuccessor(state, letter, N):
						file.write(f"  [{letter}] {successor} ")
						if letter in [2,3] or successor == 0:
							file.write("{5}\n")
						else:
							file.write("\n")
		file.write(f"--END--")





def main():
	N = int(sys.argv[1])
	assert(N > 0)
	ubaGenerator(N)


if __name__ == '__main__':
	main()