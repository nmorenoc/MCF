open (fileinput, "mcf_colloid01.dat");

$n = 0;
$t = 0.0;
$t_start = 400.0;
$t_end   = 2000.0;
$v = 0.0;
$v_tot = 0.0;

while($line = <fileinput>)
{
    @data = split(' ', $line); 

    $t = $data[1];
    $v = $data[5];
    
    if($t>=$t_start && $t<=$t_end) 
    {
	$n = $n+1;   
	$v_tot = $v_tot + $v;
    }   
}

print "v_tot     : ", $_tot, "\n";
print "n         : ", $n, "\n";
print "v_average : ", $v_tot/$n, "\n";

close(fileinput);

