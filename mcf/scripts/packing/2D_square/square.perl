$pi=3.14159265358979;
# boundary condition :
# 0: x-periodic, y-periodic/y-Lees-Edwards
# 1: x-periodic, y-wall.

$boundary = 0;

# minimum and maximum boundaryies

$min[1]    = 0.0;
$max[1]    = 32.0;
$min[2]    = 0.0;
$max[2]    = 32.0;
$VL        = ($max[1]-$min[1])*($max[2]-$min[2]);

#radius of each particle

$R        = 1.0;

# number of particles in each direction

$n[1]     = 10;
$n[2]     = 10;

$VR       = $pi * $R**2 * $n[1] * $n[2];

$con      = $VR/$VL;

# get space excluding space occupied by particles.

for ($i=1;$i<=2;$i++)
{
    $dx[$i] = ($max[$i]-$min[$i]-$n[$i]*$R);
}


# get gap between particles.
# and set initial relative position;

if( $boundary == 0) 
{
    $dx[1] = ($max[1]-$min[1]-$n[1]*2.0*$R)/$n[1];
    $dx[2] = ($max[2]-$min[2]-$n[2]*2.0*$R)/($n[2]);
    
    $x[0] = -$dx[1]/2.0 - $R;
    $y[0] = -$dx[2]/2.0 - $R;
}
else
{
    
    $dx[1] = ($max[1]-$min[1]-$n[1]*2.0*$R)/$n[1];
    $dx[2] = ($max[2]-$min[2]-$n[2]*2.0*$R)/($n[2]+1);
    
    $x[0] = -$dx[1]/2.0 - $R;
    $y[0] = -$R;
    
}

if ($dx[1] < 0.0 || $dx[2] < 0.0 ) 
{
    print "too many particles !";
    goto finished;	    
}
	

$k=0;

for ($i=1;$i<=$n[1]; $i++)
{
    $x[$i] = $x[$i-1] + $dx[1] + 2.0*$R;
    
    for ($j=1;$j<=$n[2]; $j++)
    {
	$y[$j] = $y[$j-1] + $dx[2] + 2.0*$R;
	
	$k++;
	
	$p[$k] = $x[$i];
	$q[$k] = $y[$j];
	
	print "coll_shape=1 \n";
	print "coll_radius=1,0,0.0\n";
	print 'coll_x=', $p[$k], ',', $q[$k], "\n\n";
	#print $p[$k], ' ', $q[$k], "\n";
	
    }
}


print "\nN               :", $n[1]*$n[2], "\n";
print 'volume fraction : ', $con, "\n";

 finished:
    
    
    
