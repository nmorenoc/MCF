$file_in_name = "mcf_particles00078125.out";

$num=1000;

$min_x=0.0;
$min_y=0.0;
$max_x=0.1;
$max_y=0.1;

$rc=6.0e-3;

open file_in, $file_in_name;

$i=0;
while($line=<file_in>)
    
{
    @data_in=split(' ', $line);
    $x[$i]=$data_in[0];
    $y[$i]=$data_in[1];
    $i++;    
}

$num=$i;

$k=0;
for ($i=0;$i<$num;$i++)
{
    if ($x[$i]>=$min_x+$rc &&
	$y[$i]>=$min_y+$rc &&
	$x[$i]<=$max_x-$rc &&
	$y[$i]<=$max_y-$rc)
    {
	$n[$k]=0;
	for ($j=0;$j<$num;$j++)
	{
	    if ( ($x[$j]-$x[$i])**2 + ($y[$j]-$y[$i])**2 < $rc**2)
	    {
		$n[$k]++;
	    }
	}
	$k++;
    }
}

$nk=$k;
print "num: ", $num, "\n";
print "nk : ", $nk, "\n";

for ($k=0;$k<$nk;$k++)
{
    print "n[k]: ", $n[$k], "\n";
}

close(file_in);
