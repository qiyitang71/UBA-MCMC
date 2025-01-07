dtmc

module random
	// local state
	s : [0..3] init 0;
	
	[] s=0 -> 1.0/3 : (s'=0) + 1.0/3 : (s'=1)+ 1.0/3 : (s'=2);
	[] s=1 -> 1.0/3 : (s'=0) + 1.0/3 : (s'=1)+ 1.0/3 : (s'=2);
	[] s=2 -> 1.0/3 : (s'=0) + 1.0/3 : (s'=1)+ 1.0/3 : (s'=2);
endmodule


//init s>=0 endinit