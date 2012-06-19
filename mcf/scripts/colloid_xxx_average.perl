open (fileinput, "mcf_colloid0001.dat");


$t_start = 40.0;
$t_end   = 1000.0;

$t_start = $ARGV[0];
$t_end   = $ARGV[1];

#print "t_start:", $t_start, "\n";
#print "t_end  :", $t_end, "\n";

$t = 0.0;
$n = 0;

$t_x = 0.0;
$t_x2 = 0.0;

$col  = 9;

while($line = <fileinput>)
{
    @data = split(' ', $line); 

    $t = $data[1];
    if( $t <= $t_end ) 
    {
	if( $t>=$t_start)
	{
	    $n = $n+1.0;
	    $t_x = $t_x + $data[$col];
	    $t_x2 = $t_x2 + $data[$col]**2;	
	}   
    }
    else
    {
	last;
    }
}

$t_x  = $t_x/$n;
$t_x2 = $t_x2/$n;

$t_sx = sqrt($t_x2-$t_x**2);
$t_ex = $t_sx/sqrt($n);


#print "n      : ", $n, "\n";
#print "torque, standard deviation, error\n";
#print "t_x,vari,error : ", $t_x,' ',$t_sx, ' ' , $t_ex, "\n";
print  $t_x,' ',$t_sx, ' ' , $t_ex, "\n";

close(fileinput);
