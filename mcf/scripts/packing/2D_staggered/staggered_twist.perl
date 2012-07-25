####################################################
# Staggered packing of 2D discs
# x periodic
# y wall
####################################################

$pi=3.14159265358979;

$test = $ARGV[0];

# minimum and maximum boundaryies

$min[1]    = 0.0;
$max[1]    = 16.0;
$min[2]    = 0.0;
$max[2]    = 16.0;
$VL        = ($max[1]-$min[1])*($max[2]-$min[2]);

#radius of each particle

$R        = 1.0;

# number of particles in each direction

$n[1]     = 12;
$n[2]     = 10;


# get space excluding space occupied by particles.
$shift[1]=0.2;

$shift_y[1]=0.1;
$shift_y[2]=0.05;
$shift[2]=$shift_y[1]+$shift_y[2];

#$dx[$i] = ($max[$i]-$min[$i]-2.0*$R)/$n[$i
$dx[1] = ($max[1]-$min[1]-2.0*$R-2.0*$shift[1])/($n[1]-1);
$dx[2] = ($max[2]-$min[2]-2.0*$n[2]*$R-2.0*$shift[2])/($n[2]-1);


#$x_shift[1] = $dx[1]/2.0+$R/2.0;
#$x_shift[2] = $dx[2]/2.0+$R/2.0;

$x[1] = $min[1] + $R + $shift[1];

$count = 0;

$flip = 0;

$i=1;
while($x[$i]<=$max[1]+$x_shift[1]-$R && $x[$i]<=$max[1])
{
    
    if($flip==0)
    {
	$y[1] = $min[2] + $R+ $shift[2];
	$j=1;
    
	while ($y[$j]<=$max[2]-$R)
	{
	    if ( $test == 0 ) 
	    {
		print "coll_shape=1 \n";
		print "coll_radius=1,0,0.0\n";
		print 'coll_x=', $x[$i], ',', $y[$j], "\n\n";
	    }
	    else
	    {
		print $x[$i], ' ', $y[$j], "\n";
	    }
	    $count ++;
	    $j++;
	    $y[$j] = $y[$j-1] + 2.0*$R + $dx[2];
	}
    }
    else
    {
	$y[1] = $max[2] - $R - $shift[2];
	$j=1;
	
	while ($y[$j]>=$min[2]-$R)
	{
	    if ( $test == 0 ) 
	    {
		print "coll_shape=1 \n";
		print "coll_radius=1,0,0.0\n";
		print 'coll_x=', $x[$i], ',', $y[$j], "\n\n";
	    }
	    else
	    {
		print $x[$i], ' ', $y[$j], "\n";
	    }
	    $count ++;
	    $j++;
	    $y[$j] = $y[$j-1] - 2.0*$R - $dx[2];
	}
    }
    
    $i++;
    $x[$i] = $x[$i-1] + $dx[1] + 2.0*$R;
    
    if ($flip ==0 ) 
    { 
	$flip = 1;
    }
    else
    {
	$flip =  0;
    }
    
}

$VR       = $pi * $R**2 * $count;
$con      = $VR/$VL;

if ($test == 0 )
{
    print "\nN               :",$count, "\n";
    print 'volume fraction : ', $con, "\n";
}

#print $dx[1], ' ', $dx[2], "\n";	

 finished:
    
    
    
