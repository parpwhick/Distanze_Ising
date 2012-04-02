#!/usr/bin/perl


for($beta=0.39;  $beta<=0.46;$beta=$beta+0.01){
	$beta=sprintf("%.2f",$beta);
	open dove, ">>beta_$beta";
	$oldval=0;
	for($l=10;$l<400;$l=$l * 1.2){
#	for($l=1;$l<10;$l=$l +1){
		$l=sprintf("%d",$l);
		open risultato, "./ising -seqnum 20 -lattice -lato $l -beta $beta|";
		$valori=<risultato>;
		@a=split " ",$valori;
		print  "b: $beta, $valori";
		print dove $valori;
		last if ($l> 200 && $a[1]>5 && $a[1]-$oldval < 0.005 && $a[2] < 0.1);
		$oldval=$a[1];
	}
}
