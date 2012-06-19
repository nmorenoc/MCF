
$a  = 0.02;
$n  = 7;
$l0 = 0.06;
$dl = 0.01;
$pi = 3.14159265358979;
$U  = 1.78550706649452e-4;
$ksai = 1.0e-1;
$rho = 1.0e+3;
$Csound = 1.0e-2;

print 'l c ,Q, U, F, f, Re, Csound, Mach', "\n";

for ($i=0;$i<7;++$i)
{
    $l[$i] = $l0+$i*$dl;
    print $l[$i], ' ';    
    $c[$i] = 4.0*$pi*$a**3 / 3.0/$l[$i]**3;
    print $c[$i],' ';
    $Q[$i] = 1 - 1.7601 * $c[$i]**(1.0/3.0) + $c[$i] - 1.5593 * $c[$i] **2;
    print $Q[$i], '  ';
    print $U, ' ';
    $F[$i] = 6.0*$pi*$a*$ksai*$U/$Q[$i];
    print $F[$i], ' ' ;
    $f[$i] = $F[$i]/$rho/($l[$i]**3 -4.0*$pi*$a**3 / 3.0);
    print $f[$i], ' ';
    print $U*$a/($ksai/$rho), ' ';
    print $Csound, ' ';
    print $U/$Csound, ' ';
    print "\n";
}


