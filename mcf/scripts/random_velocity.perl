
$pi = 3.14159265358979;
$nx = 32;
$ny = 32;

$v0 = 1.0;

$lx = 1.0;
$ly = 1.0;


$hx = $lx/$nx;
$hy = $ly/$ny;

$rho = 1.0;

$m  = $rho * $hx * $hy;

$file_name_a = ">random_a.dat";
$file_name_b = ">random_b.dat";
$file_name_c = ">random_c.dat";


$Mom[1] = 0.0;
$Mom[2] = 0.0;

my @vx;
my @vy;

#print @vx;
#print @vy;

open file_a, $file_name_a;
open file_b, $file_name_b;
open file_c, $file_name_c;

print file_a $nx*$ny, "\n";
print file_b $nx*$ny, "\n";

$id = 0;

for ($j=1; $j<=$ny;++$j)
{    
    for ($i=1; $i<=$nx;++$i)
    {
	$myrand = rand();
	$vx[$i-1 + $nx*($j-1)] = $v0*cos(2*$pi*$myrand);
	$vy[$i-1 + $nx*($j-1)] = $v0*sin(2*$pi*$myrand); 
	#$vx[$i-1 + $nx*($j-1)] = 0.0;
	#$vy[$i-1 + $nx*($j-1)] = 0.0; 

	$Mom[1] = $Mom[1] + $vx[$i-1+$nx*($j-1)];
	$Mom[2] = $Mom[2] + $vy[$i-1+$nx*($j-1)];

    }
}
    

$Mom[1] = $Mom[1] / ($nx * $ny);
$Mom[2] = $Mom[2] / ($nx * $ny);

#print $Mom[1], ' ', $Mom[2], "\n";

for ($j=1; $j<=$ny;++$j)
{    
   $y = ($j-0.5) * $hy;
   
   for ($i=1; $i<=$nx;++$i)
   {
       $x = ($i-0.5) * $hx;

       $id = $id + 1;
       print file_a  $x, ' ', $y, ' ', $vx[$i-1+$nx*($j-1)], ' ', $vy[$i-1+$nx*($j-1)], ' ',$rho, ' ', $m, ' ', $id, ' ', 0.0, "\n";
       print file_b  $x, ' ', $y, ' ', $vx[$i-1+$nx*($j-1)]-$Mom[1], ' ', $vy[$i-1+$nx*($j-1)]-$Mom[2], ' ',$rho, ' ', $m, ' ', $id, ' ', 0.0, "\n";
       print file_c  $x, ' ', $y, ' ', $vx[$i-1+$nx*($j-1)]-$Mom[1], ' ', $vy[$i-1+$nx*($j-1)]-$Mom[2], ' ',$rho, ' ', $m,  "\n";
   }
}

close(file_a);
close(file_b);
close(file_c);






