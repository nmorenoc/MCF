open (fileinput, "colloids_sid.dat");
open (fileoutput,">mcf_colloid0001.dat");

$step_freq = 500;
$dt        = 0.461731678571429;

$step_freq= $ARGV[0];
$dt       = $ARGV[1];


print "step_freq:", $step_freq, "\n";
print "dt       :", $dt, "\n";


$n=0;
while($line = <fileinput>)
{
    @data = split(' ', $line); 
    
    $step=$n*$step_freq;
    $time=$n*$dt;

    print fileoutput  $step, ' ', $time;
    for ($i=0;$i<=4;$i++)
    {
	print fileoutput  ' ', $data[$i];
    }
    print fileoutput "\n";
    $n++;
}


close(fileinput);
close(fileoutput);
