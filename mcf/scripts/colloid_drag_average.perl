open (fileinput, "mcf_colloid0001.dat");

$n = 0;
$t = 0.0;
$t_start = 200.0;
$t_end   = 1000.0;

$t_start = $ARGV[0];
$t_end   = $ARGV[1];

#print "t_start:", $t_start, "\n";
#print "t_end  :", $t_end, "\n";

$drag = 0.0;
$drag_tot = 0.0;


while($line = <fileinput>)
{
    @data = split(' ', $line); 

    $t = $data[1];
    $drag = $data[2];
    
    if ($t<=$t_end)
    {
	if($t>=$t_start)
	{
	    $n = $n+1;   
	    $drag_tot = $drag_tot + $drag;
	}   
    }
    else
    {
	last;
    }

    #print $t, ' ', $drag, ' ', "\n";
}

#print "drag_tot     : ", $drag_tot, "\n";
#print "n            : ", $n, "\n";
#print "drag_average : ", $drag_tot/$n, "\n";
print $drag_tot/$n, "\n";

close(fileinput);

