#!/usr/bin/perl


for($beta=100000; $beta<=100000;$beta=$beta+10000){
	$beta=sprintf("%.2f",$beta);
	open dove, ">>beta_$beta";
	$oldval=0;
	
	for($l=10;$l<10;$l=$l +1){
		$l=sprintf("%d",$l);
		open risultato, "./ising -seqnum 100 -seqlength $l -beta $beta|";
		$valori=<risultato>;
		@a=split " ",$valori;
		print  "$beta: $valori";
		print dove $valori;
	}
	
	$oldval=0;
	#for($l=10;$l<300000;$l=$l * 1.3){
	for($l=300000;$l<=300000;$l=$l * 1.3){
		$l=sprintf("%d",$l);
		open risultato, "./ising -seqnum 30 -seqlength $l -beta $beta|";
		$valori=<risultato>;
		@a=split " ",$valori;
		print  "$beta: $valori";
		print dove $valori;
		last if ($a[1]>1.1); # && $a[1]-$oldval < 0.005);
		$oldval=$a[1];
	}

}
