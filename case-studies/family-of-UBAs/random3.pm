dtmc

module random
	// local state
	s : [0..999] init 0;
	
	[] s<997 -> 0.25 : (s'=s) + 0.25 : (s'=s+1)+ 0.25 : (s'=s+2)+ 0.25 : (s'=s+3);
	[] s=997 -> 0.25 : (s'=s) + 0.25 : (s'=s+1)+ 0.25 : (s'=s+1)+ 0.25 : (s'=0);
	[] s=998 -> 0.25 : (s'=0) + 0.25 : (s'=1)+ 0.25 : (s'=s+1)+ 0.25 : (s'=s);
	[] s=999 -> 0.25 : (s'=0) + 0.25 : (s'=1)+ 0.25 : (s'=2)+ 0.25 : (s'=s);
endmodule


//init s>=0 endinit
label "sigma_0" = s=0;
label "pi_0" = s=1;
label "hash_0" = s=2;
label "dollar_0" = s=3;
