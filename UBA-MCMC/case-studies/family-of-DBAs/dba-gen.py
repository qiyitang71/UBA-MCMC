#!/usr/bin/env python3
import sys
"""
Alphabet = {0,1,2} 
0: sigma
1 : pi
2 : hash
3 : dollar

N is the no of boxes.

Each state of the automata is a number in range (0,2^N -2). The state 0 represents the state q_$. For states > 1
the binary representation of this number is a binary string of length n, and the value of bit in i-th position stores if q_i is in the state. 
A value of 1 denotes that q_i is in the state. For eg : for n = 3, 011
denotes the state {q_0,q_1}. 

Accepting transition are all labelled with 0. 
"""

def findSuccessor(state, letter, N): #gives list of successor states on reading letter from state
	if letter == 0: # sigma
		return([((((state+1) >> (N-1))&1)|(2**N-1)&((state+1) << 1))-1])

	elif letter == 1: #pi
		lastTwoBit = (state + 1) & 3
		if lastTwoBit in [0,3]:
			return([state])
		else:
			return([((state+1) ^ 3) - 1])

	elif letter ==2: # hash
		if state == 0:
			return([2**N - 3 ])
		else:
			return([(state+1) - ((state+1) & 1) - 1])

def labelGenerator(letter):
	labelValues = []
	for lt in range(3):
		if lt != letter:
			labelValues.append("!"+str(lt))
		else:
			labelValues.append(str(lt))
	label = " & ".join(labelValues)
	return(label)

def countBit(n: int) -> int:
    return sum(b == '1' for b in bin(n))
    
def ubaGenerator(N):
	with open(f"DBAs/dba-{N}.hoa",'w') as file:
		file.write(f"HOA: v1\n")
		file.write(f"name: \"dba-{N}\"\n")
		NO_OF_STATES = 2**N -2
		file.write(f"States: {NO_OF_STATES}\n")
		file.write(f"Start: 0\n")
		file.write(f"acc-name: Buchi\n")
		file.write(f"Acceptance: 1 Inf(0)\n") #each acc transition is given label 0
		file.write(f"AP: 3 \"sigma\" \"pi\" \"hash\"\n")
		file.write("properties: trans-labels explicit-labels trans-acc unambiguous\n")
		file.write(f"--BODY--\n")
		for state in range(NO_OF_STATES):
			file.write(f"State: {state}\n")
			for letter in range(3):
				for successor in findSuccessor(state, letter, N):
					file.write(f"  [{labelGenerator(letter)}] {successor} ")
					if (countBit(state+1) != countBit(successor+1) and letter == 2):
						file.write("{0}\n")
					else:
						file.write("\n")
		file.write(f"--END--")



def main():
	N = int(sys.argv[1])
	assert(N > 0)
	ubaGenerator(N)


if __name__ == '__main__':
	main()