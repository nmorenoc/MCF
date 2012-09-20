###################################################
#Calculate average of a set of data
###################################################

$file_name=$ARGV[0];
open (file_in, $file_name);


###################################################
# time period of interest.
###################################################
$t_start = 2.0e2;
$t_end   = 1.0e6;

$t_start = $ARGV[1];
$t_end   = $ARGV[2];

#print "t_start:", $t_start, "\n";
#print "t_end  :", $t_end, "\n";

$t = 0.0;
$n = 0;

###################################################
# number of columns of data NC,
# interest columns [N1,N2]
###################################################
$N1=$ARGV[3];
$N2=$ARGV[4];

###################################################
# number of rows interested
###################################################
$NR=0;
for($i=0;$i<=$N2-$N1;$i++)
{
    $average[$i]=0.0;
}

while($line = <file_in>)
{
    @data = split(' ', $line); 
    
    $step=$data[0];
    $time=$data[1];
    
    if ($time<$t_start)
    {
	next;
    }
    elsif($time<=$t_end)
    {
	for($i=$N1;$i<=$N2;$i++)
	{
	    $average[$i-$N1]+=$data[$i];
	}
	$NR++;
    }
    elsif ($time>$t_end)
    {
	last;
    }
}

for ($i=0;$i<=$N2-$N1;$i++)
{
    $average[$i]/=$NR;
    print $average[$i], " ";
}
print "\n";

close(file_in);
    
    
