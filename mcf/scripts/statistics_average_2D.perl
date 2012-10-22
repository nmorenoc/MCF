###################################################
# Calculate average statistics in 2D
###################################################

open (file_in, "mcf_statistic.dat");

#########################################
# time period of interest.
#########################################
$t_start = 40.0;
$t_end   = 1000.0;

$t_start = $ARGV[0];
$t_end   = $ARGV[1];

#print "t_start:", $t_start, "\n";
#print "t_end  :", $t_end, "\n";


$t = 0.0;
$n = 0;

#########################################
#columns of interest
#########################################
$c_start=3;
$c_end=3;

$c_start = $ARGV[2];
$c_end   = $ARGV[3];
$num_c=$c_end-$c_start+1;

for ($i=0;$i<$num_c;$i++)
{
    $aver[$i] =0.0;
    $aver2[$i]=0.0;
}

while($line = <file_in>)
{
    @data = split(' ', $line); 
    
    $t = $data[1];
    if( $t<=$t_end )
    {
	if( $t>=$t_start )
	{
	    $n = $n+1.0;
	    for ($i=$c_start;$i<=$c_end;$i++)
	    {
		$aver[$i-$c_start] += $data[$i]; 
		$aver2[$i-$c_start] += $data[$i]**2; 
	    }	    
	}   
    }
    else
    {
	last;
    }
}

########################################
#n=0: indicate no data for processing.
########################################
if ( $n==0 )
{
    print "n=0! \n";
}
else
{
    for ($i=0;$i<$num_c;$i++)
    {
	$aver[$i] /= $n;
	$aver2[$i] /= $n;
	$avers[$i]=sqrt($aver2[$i]-$aver[$i]**2);
	$avere[$i]=$avers[$i]/sqrt($n);
	print "", $aver[$i], ' ', $avers[$i], ' ', $avere[$i], ' ', $n, "\n";
    }   
    
} # n=0

close(file_in);
    
    
