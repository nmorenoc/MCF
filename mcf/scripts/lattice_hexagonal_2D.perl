open file_out, ">hexagonal_lattice.dat";

$h = 0.002;
$dx = $h* sqrt(3.0)/2.0;
$L  = 0.1;
$prec = 1.0e-10;
$rho = 1000.0;
$m_tot = $rho * $L**2;
$n = 0;
$m = 0;
$nx = 0;
$ny = 0;

$x = $dx/2.0;
$x_coef = 0;

while ($x<$L) 
{

    ++$nx ;
    $y = (2-$x_coef)* $h/2.0;
    $y_coef = $x_coef;
    
    $ny = 0;
    while($y<$L)
    {
	++$ny;
	++$n;
	$y_coef = ($y_coef + 1) % 2;
	#$x = $x + $h* (2-$x_coef);
	$y= $y+$h;
    }
   # print "nx : ", $nx, "\n";
    $x = $x + $dx;
    $x_coef = ($x_coef + 1) % 2;
}

print "nx : ", $nx, "\n";
print "ny : ", $ny, "\n";


$m = $m_tot / $n;
$x = $dx/2.0;
$x_coef = 0;

print file_out $n, "\n";

while ($x<$L) 
{

    $y = (2-$x_coef)* $h/2.0;
    $y_coef = $x_coef;
    
    while($y<$L)
    {
	#print $x,' ', $y, ' ', 0.0, ' ', 0.0, ' ', 0.0, ' ', $m, ' ', 0,' ', 0, "\n";
	print file_out $x,' ', $y, ' ', 0.0, ' ', 0.0, ' ', 0.0, ' ', $m, ' ', 0,' ', 0, "\n";	

	$y_coef = ($y_coef + 1) % 2;
	#$x = $x + $h* (2-$x_coef);
	$y= $y+$h;

    }
    $x = $x + $dx;
    $x_coef = ($x_coef + 1) % 2;

}

print "num_part : ", $n, "\n";

close(file_out);
