
#----------------------------------
# Generate a file which contains
# theoretical solution of couette
# flow, Morris J.P. et al. 1997.
# The discrete points position is
# on simple evenly distributed grid.
#----------------------------------

$epsilon = 1.0e-20;
$pi      = 3.14159265358979;

$npart = 100;
$t   = 0.01;
$mu  = 8.46;
$L   = 2.0;
$rho = 1.0;
$v0  = 0.423;
$dh   = $L/$npart;

$timestring = substr($t, 2);
$filename = "couette" . '_' . $npart . '_'.  $timestring .  ".dat";
#$filename = "couette" . '_' . $npart . '_'.  20 .  ".dat";

open file,'>',$filename or die ;


for ($y = $dh/2.0; $y < $L; $y=$y+$dh)
{
    $Vx = $v0*$y/$L ;
    
    $n = 1;
    $temp=0.0;
    
    do
    {
	$temp = 2.0 * $v0/$n/$pi;
	$temp = $temp * sin($n * $pi / $L * $y);
	$temp = $temp * exp(-$mu * ($n**2) * ($pi**2) / ($L**2)  *$t );
	
	if (abs($temp) > $epsilon)
	{
	    if($n % 2 == 0)
	    {
		$Vx = $Vx + $temp;
	    }
	    else
	    {
		$Vx = $Vx - $temp;
	    }
	}
	
	$n = $n +1;
    }
    while(abs($temp) > $epsilon);

    print file $y,' ', $Vx, "\n";
}

close(file);

