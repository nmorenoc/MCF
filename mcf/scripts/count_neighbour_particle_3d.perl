$file_in_name = "mcf_particles00000000.out";

$num=1000;

$min_x=0.0;
$min_y=0.0;
$min_z=0.0;

$max_x=0.2;
$max_y=0.2;
$max_z=0.2;

$rc=6.0e-3;

open file_in, $file_in_name;

$i=0;

while($line=<file_in>)
    
{
    @data_in=split(' ', $line);
    $x[$i]=$data_in[0];
    $y[$i]=$data_in[1];
    $z[$i]=$data_in[2];
    $i++;    
}

$num=$i;
print "num: ", $num, "\n";

$k=0;
for ($i=0;$i<$num;$i++)
{
    if ($x[$i]>=$min_x+$rc &&
	$y[$i]>=$min_y+$rc &&
	$z[$i]>=$min_z+$rc &&
	$x[$i]<=$max_x-$rc &&
	$y[$i]<=$max_y-$rc &&
	$z[$i]<=$max_z-$rc )
    {
	$n[$k]=0;
	for ($j=0;$j<$num;$j++)
	{
	    if ( ($x[$j]-$x[$i])**2 + ($y[$j]-$y[$i])**2 
		 + ($z[$j]-$z[$i])**2 < $rc**2)
	    {
		$n[$k]++;
	    }
	}
	$k++;
	if($k>10)
	{
	    last;
	}
    }
}

$nk=$k;
print "nk : ", $nk, "\n";
$n_aver=0;

for ($k=0;$k<$nk;$k++)
{
    print "n[k]: ", $n[$k], "\n";
    $n_aver+=$n[$k]/$nk;
}

print "n_average: ", $n_average, "\n";

close(file_in);
