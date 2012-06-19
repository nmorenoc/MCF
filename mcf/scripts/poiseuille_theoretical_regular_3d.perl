
$epsilon = 1.0e-18;
$pi      = 3.14159265359;

$npart = 10;
$t   = 0.1;
$mu  = 1.0e-6;
$L   = 1.0e-3;
$rho = 1.0e+3;
$F   = 1.0e-4;
$dh   = $L /$npart;
$temp = 0;
$n    = 0;
$timestring = substr($t, 2);
$filename = "poiseuille" . '_' . $L/$dh . '_'.  $timestring .  ".dat";
open file,'>',$filename or die ;


for ($y = $dh/2.0; $y < $L; $y=$y+$dh)
{
    
    $Vx = $F/(2.0*$mu) * $y *($y-$L);
    
    $n = 0;
    $temp=0.0;
    do
    {
	$temp = 4.0*$F*($L**2)/( ($mu*($pi**3)) * (2*$n + 1)**3 );
	$temp = $temp * sin($pi * $y / $L * (2*$n+1));
	$temp = $temp * exp( -(2*$n+1)**2*($pi**2)*$mu *$t/($L**2) );
	if ($temp > $epsilon)
	{
	    $Vx = $Vx + $temp;
	}
	$n = $n +1;
    }
    while($temp > $epsilon);

    print file $y,' ', -$Vx, "\n";
}

    close(file);

