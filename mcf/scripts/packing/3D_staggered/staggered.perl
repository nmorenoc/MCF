####################################################
# Staggered packing of 3D spheres.
# x, y periodic and z solid wall.
####################################################

$pi=3.14159265358979;

$test = $ARGV[0];

# dimension of the problem.

$dim  = 3;

# minimum and maximum boundaryies.

$min[1]    = 0.0;
$max[1]    = 7.0;
$min[2]    = 0.0;
$max[2]    = 7.0;
$min[3]    = 0.0;
$max[3]    = 7.0;

$VL        = ($max[1]-$min[1])*($max[2]-$min[2])*($max[3]-$min[3]);

# radius of each particle.

$R        = 1.0;

# number of locations for particles in each direction.

$n[1]     = 3;
$n[2]     = 3;
$n[3]     = 6;

# get space for particles away from boundaries.

$shift[1]=0.1;
$shift[2]=0.1;
$shift[3]=0.1;

# calculate space between neigboring particles in each direction.

# for periodic boundary and solid wall boundary.
for ($i=1;$i<=$dim;$i++)
{
    $dx[$i] = ($max[$i]-$min[$i]-2.0*$shift[$i])/$n[$i];
}

# for solid wall boundary.
#$dx[3] = ($max[3]-$min[3]-2.0*$R-2.0*$shift[3])/($n[3]-1);


print "dx: ", $dx[1], ' ', $dx[2], ' ', $dx[3], "\n\n";


# reset counter of particles.

$count = 0;

# first location of particle center in y direction

$y[1]  = $min[2] + $dx[2]/2.0 + $shift[2];
$j     = 1;

#  particle's y location must be smaller than max[2] and
#  particle's sufrace must be smaller than another particle
#  at periodic boundary with a shift[1].

while( $y[$j]<=$max[2] && $y[$j]+$R <= $max[2]+$shift[2] )
{
    # first location of particle center in x direction.
    
    $x[1] = $min[1]  + $dx[1]/2.0 + $shift[1];
    $i    = 1;
    
    # flipping indicator.

    $flip = 0;
    
   #  particle's x location must be smaller than max[1] and
   #  particle's sufrace must be smaller than another particle
   #  at periodic boundary with a shift[1].
    
    while($x[$i]<=$max[1] && $x[$i]+$R <= $max[1]+$shift[1] )
    {
	# first location of particle center in z direction.
	
	$z[1] = $min[3] + $shift[3]+ $dx[3] + $flip * $dx[3];
	
	$k=1;
	
	#  particle's surface must be smaller than max[3]
	
	while ($z[$k]<=$max[3]-$R)
	{
	    if ( $test == 0 ) 
	    {
		print "coll_shape=2 \n";
		print "coll_radius=1,0,0.0,0.0\n";
		print 'coll_x=', $x[$i], ',', $y[$j], ',', $z[$k],"\n\n";
	    }
	    else
	    {
		print $x[$i], ' ', $y[$j], ' ', $z[$k],"\n";
	    }
	    
	    $count ++;
	    $k++;
	    $z[$k] = $z[$k-1] + 2*$dx[3];

	} # z
	
	$i++;
	$x[$i] = $x[$i-1] + $dx[1];
	
	if ( $flip ==0 ) 
	{ 
	    $flip = 1;
	}
	else
	{
	    $flip =  0;
	}
	
    } # x
    $j++;
    $y[$j] = $y[$j-1] + $dx[2];
    
} # y

$VR   = $count* 4.0*$pi*$R**3/3.0 ;
$con  = $VR/$VL;

if ($test == 0 )
{
    print "\nN               :",$count, "\n";
    print "VR                :", $VR, "\n";
    print "VL:               :", $VL, "\n";
    print "volume fraction   :", $con, "\n";
}


 finished:
    
    
    
