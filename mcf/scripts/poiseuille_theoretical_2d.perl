open file_n, "thin.dat";
open file_a, ">poiseuille_b.dat";

$t   = 10;
$epsilon = 1.0e-18;
$pi      = 3.14159265359;
$mu  = 0.06;
$L   = 1.0;
$rho = 1.0;
$F   = 0.0135765;
$temp = 0;
$n    = 100;

while($line_n=<file_n>)
{
    
    @data_n = split(' ', $line_n);
    
    $y = $data_n[1];
    
    $Vx = $F/(2.0*$mu) * $y *($y-$L);
    
    $n = 0;
    $temp=0.0;
    do
    {
	$temp = 4.0*$F*($L**2)/( ($mu*($pi**3)) * (2*$n + 1)**3 );
	$temp = $temp * sin($pi * $y / $L * (2*$n+1));
	$temp = $temp * exp( -(2*$n+1)**2*($pi**2)*$mu *$t/($L**2) );
	
	$Vx = $Vx + $temp;	
	$n = $n +1;
    }
    while($temp > $epsilon);

    #print $y, ' ', -$Vx, "\n";
    print file_a  $y,' ', -$Vx, "\n";
}

close(file_n);
close(file_a);


