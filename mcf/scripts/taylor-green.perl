
$step_end  = 0;
$dt = 0.1;
$pi = 3.14159265358979;

$nx = 32;
$ny = 32;

$mu = 0.1;


$a  = 0.0016;
$lx = 2*$pi;
$ly = 2*$pi;

$Re = 1.0;

$k  = 1;

$hx = $lx/$nx;
$hy = $ly/$ny;

$rho = 1.0;

$m  = $rho * $hx * $hy;

$file_name = ">taylor-green";
$file_name = $file_name . '_' . $step_end . ".dat";

open file_a, $file_name;

print file_a $nx*$ny, "\n";

for ($step = 0; $step<=$step_end; ++$step)
{
    $t = $step * $dt;
    
    for ($j=1; $j<=$ny;++$j)
    {    
	$y = ($j-0.5) * $hy;
	
	for ($i=1; $i<=$nx;++$i)
	{
	    $x = ($i-0.5) * $hx;
	    $vx = - sqrt($a)* cos($k*$x)*sin($k*$y) * exp(-(2*$k)**2*$t/$Re);
	    $vy =  sqrt($a) * sin($k*$x)*cos($k*$y) * exp(-(2*$k)**2*$t/$Re);
	    #print $x, ' ', $y, ' ', $vx, ' ', $vy, ' ',$rho, ' ', $m, ' ', 0.0, ' ', 0.0, "\n";
	    print file_a  $x, ' ', $y, ' ', $vx, ' ', $vy, ' ',$rho, ' ', $m, ' ', 0.0, ' ', 0.0, "\n";
	}
    }
    
}
close(file_a)






