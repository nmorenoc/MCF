open file_n, "thin.dat";
open file_a, ">poiseuille_b.dat";

$t   = 0.1;
$epsilon = 1.0e-18;
$pi      = 3.14159265359;
$mu  = 1.0e-6;
$L   = 1.0e-3;
$rho = 1.0e+3;
$F   = 1.0e-4;
$temp = 0;
$n    = 0;

while($line_n=<file_n>)
{
    
    @data_n = split(' ', $line_n);
    
    $y = $data_n[2];
    
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


