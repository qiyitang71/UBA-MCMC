dtmc

module random
	// local state
	s : [0..19] init 0;
	
	[] s<17 -> 0.25 : (s'=s) + 0.25 : (s'=s+1)+ 0.25 : (s'=s+2)+ 0.25 : (s'=s+3);
	[] s=17 -> 0.25 : (s'=s) + 0.25 : (s'=s+1)+ 0.25 : (s'=s+2)+ 0.25 : (s'=0);
	[] s=18 -> 0.25 : (s'=0) + 0.25 : (s'=1)+ 0.25 : (s'=s+1)+ 0.25 : (s'=s);
	[] s=19 -> 0.25 : (s'=0) + 0.25 : (s'=1)+ 0.25 : (s'=2)+ 0.25 : (s'=s);
endmodule


//init s>=0 endinit
label "sigma" = mod(s,4)=0;
label "pi" = mod(s,4)=1;
label "hash" = mod(s,4)=2;
label "dollar" = mod(s,4)=3;