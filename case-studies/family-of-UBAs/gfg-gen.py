#!/usr/bin/env python3
import sys
"""
Alphabet = {0,..,7} 
0 : <sigma,0>
1 : <pi,0>
2 : <hash,0>
3 : <dollar,0>
4 : <sigma,1>
5 : <pi,1>
6 : <hash,1>
7 : <dollar,1>

N is the no of boxes.

Each state of the automata is a number in range 0 to N+1. The state N represents the state q$ and the state N+1 
representes the sink state(bottom in paper). A state in i range [0,N-1] denotes the state q_i. 

Accepting transition are all labelled with 0. 
"""

def findSuccessor(state, letter, N): #gives list of successor states on reading letter from state
	if state == N+1: #sink state
		return([N+1])
	elif state == N: # q$
		if letter%4 < 3: # any letter other than dollar goes to sink
			return([N+1])
		elif letter == 3:
			return([N])
		else:
			assert(letter == 7)
			return([0])
	elif letter%4 == 3: # from any q_i reading dollar goes to sink
		assert(state < N)
		return([N+1])
	elif letter//4 ==0: # <a,0> letters from q_1 goes to dollar
		assert(state < N)
		return([N])

	elif letter == 4: # <sigma,1> 
		return([(state+1)%N])

	elif letter == 5: #<pi,1>
		if state > 1:
			return([state])
		else:
			return([(state+1)%2])
	else:
		assert(letter == 6) # <hash,1>
		if state == 0:
			return(list(range(N)))
		else:
			return([state])



def labelGenerator(letter):
	labelValues = []
	for lt in range(8):
		if lt != letter:
			labelValues.append("!"+str(lt))
		else:
			labelValues.append(str(lt))
	label = " & ".join(labelValues)
	return(label)

def gfgGenerator(N):
	with open(f"GFGs/gfg-{N}.hoa",'w') as file:
		file.write(f"HOA: v1\n")
		file.write(f"name: \"gfg-{N}\"\n")
		NO_OF_STATES = N+2
		file.write(f"States: {NO_OF_STATES}\n")
		file.write(f"Start: 0\n")
		file.write(f"acc-name: Buchi\n")
		file.write(f"Acceptance: 1 Inf(0)\n") #each acc transition is given label 0
		file.write("AP: 8 \"sigma_0\" \"pi_0\" \"hash_0\" \"dollar_0\" \"sigma_1\" \"pi_1\" \"hash_1\" \"dollar_1\"\n")
		file.write("properties: trans-labels explicit-labels trans-acc\n")
		file.write(f"--BODY--\n")
		for state in range(NO_OF_STATES):
			file.write(f"State: {state}\n")
			for letter in range(8):
				for successor in findSuccessor(state, letter, N):
					file.write(f"  [{labelGenerator(letter)}] {successor} ")
					if (state,letter) in [(0,6), (N,3),(N,7)]:
						file.write("{0}\n")
					elif successor == N:
						file.write("{0}\n")
					else:
						file.write("\n")
		file.write(f"--END--")





def main():
	N = int(sys.argv[1])
	assert(N > 0)
	gfgGenerator(N)


if __name__ == '__main__':
	main()