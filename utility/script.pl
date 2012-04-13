#!/usr/bin/perl


for($beta=1; $beta<=1;$beta=$beta+10){
	$beta=sprintf("%.2f",$beta);
	open dove, ">>beta_$beta";
	$oldval=0;
	
	for($l=1;$l<10;$l=$l +1){
		$l=sprintf("%d",$l);
		open risultato, "./ising -seqnum 100 -seqlength $l -beta $beta|";
		$valori=<risultato>;
		@a=split " ",$valori;
		print  "$beta: $valori";
		print dove $valori;
	}
	
	$oldval=0;
	for($l=10;$l<30000;$l=$l * 1.1){
		$l=sprintf("%d",$l);
		open risultato, "./ising -seqnum 100 -seqlength $l -beta $beta|";
		$valori=<risultato>;
		@a=split " ",$valori;
		print  "$beta: $valori";
		print dove $valori;
		last if ($a[1]>1.28); # && $a[1]-$oldval < 0.005);
		$oldval=$a[1];
	}

}
