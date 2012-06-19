
$n = 0;
$t = 0.0;

$file_name='';
$x_start = 0.0;
$x_end   = 0.0;

$file_name = $ARGV[0];
$x_start   = $ARGV[1];
$x_end     = $ARGV[2];

#$print "file_name: ", $file_name, "\n";
#print "x_start:", $x_start, "\n";
#print "x_end  :", $x_end, "\n";

open (fileinput, $file_name);

$v = 0.0;
$v_tot = 0.0;


while($line = <fileinput>)
{
    @data = split(' ', $line); 

    $x = $data[0];
    $v = $data[2];
    
    if ($x<=$x_end)
    {
	if($x>=$x_start)
	{
	    $n = $n+1;   
	    $v_tot = $v_tot + $v;
	}   
    }   

    #print $x, ' ', $v, ' ', "\n";
}

#print "v_tot     : ", $v_tot, "\n";
print "n            : ", $n, "\n";
#print "v_average : ", $v_tot/$n, "\n";
print $v_tot/$n, "\n";

close(fileinput);

