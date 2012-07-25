####################################################
# read in colloid file from previous run
####################################################

$pi=3.14159265358979;

$test = $ARGV[0];
$name_in=$ARGV[1];

open file_in, $name_in;

#radius of each particle

$R        = 1.0;

# number of particles in each direction

$count=0;

while ( $line=<file_in>)
{
    @data=split(' ',$line);

    $x=$data[0];
    $y=$data[1];
    $sid=$data[7];

    if ( $sid == 0 ) 
    {
	if ( $test == 0 ) 
	{
	    print "coll_shape=1 \n";
	    print "coll_radius=1,0,0.0\n";
	    print 'coll_x=', $x, ',', $y, "\n\n";
	}
	else
	{
	    print $x, ' ', $y, "\n";
	}
	$count ++;
    }
}

print "count: ", $count, "\n";
   
    
close(file_in);
