use Math::Trig;

$file_in="mcf_colloid0001.dat";
$file_out=">mcf_colloid.dat";

open (fileinput, $file_in);
open (fileoutput, $file_out);


$t_start = 0.0;
$t_end   = 1000.0;

#$t_start = $ARGV[0];
#$t_end   = $ARGV[1];

print "Processing ", $file_in, "...\n";
#print "t_start : ", $t_start, "\n";
#print "t_end   : ", $t_end, "\n";

$col= 5;
$s  = 0;
$t  = 0.0;
$temp =0.0;
$cos_th_old = 0.0;
$cos_th = 0.0;
$sin_th = 0.0;
$tan_th = 0.0;
$th =0.0;
$omega=0.0;

$n = 0;

$cos_th_old = 1.0;

while($line = <fileinput>)
{
    @data = split(' ', $line); 

    $t = $data[1];
    
    if( $t<=$t_end) 
    {
	if($t>=$t_start) 
	{
	    $s = $data[0];
	    
	    $n = $n+1;
	    
	    # sum of the trace of rotation matrix.
	    
	    $temp=$data[$col];
	    
	    for ($i=4;$i<=8;$i+=4)
	    {
		$temp += $data[$col+$i];
	    }
	    
	    # compute rotation angle using the trace.
	    $cos_th = ($temp-1.0)/2.0;
	    # if cos_th starts to increase, it means
            # we are at 3rd, 4th quadrant and 
	    # sin should be negative. 
	    if ($cos_th > $cos_th_old)
	    {
		$sin_th = -sqrt(abs(1-$cos_th**2));
	    }
	    else
	    {
		$sin_th = sqrt(abs(1-$cos_th**2));
	    }
            # use atan to compare with Jeffery 1922 solution
            # easily.
	    $tan_th = $sin_th/$cos_th;
	    $th = atan($tan_th);
	    $cos_th_old = $cos_th;
            # this is the angular speed rotating around y-aixs
	    $omega = $data[15];

	    print  fileoutput $s,' ',$t, ' ', $th, ' ', $omega, "\n";
	}
    }
    else
    {
	last;
    }
}

print $n, " time steps of data! \n";
close(fileinput);
close(fileoutput);
