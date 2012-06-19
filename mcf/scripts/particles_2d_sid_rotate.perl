#print $ARGV[0];
#print $ARGV[1];
###################################################
# Mark species ID plus one extra marker to
# indicate rotation of a colloidal particle.
###################################################
$ARGV[1] = '>>' . $ARGV[1];
open file_in, $ARGV[0];
open file_out, $ARGV[1];

#suppose maximum we have 200 colloid
#suppose maximum 1 billion particles

$num_colloid=0;
for ($i=1;$i<=200;$i++)
{
    $pid[$i]=1.e9;
    $row[$i]=0;
}

$col=10;
$n=0;

while($line = <file_in>)
{
 
    @data_in = split(' ', $line);
    
    #check for each colloidal boundary particle
    if ($data_in[7] > 0)
    {
# check the biggest particle id of each
# colloidal boundary particle in this colloid
# and record its row number.
	if ($data_in[6]<$pid[$data_in[7]])
	{
	    $pid[$data_in[7]]=$data_in[6];
	    $row[$data_in[7]]=$n;
	}
# check the maximum index of colloid.
	if ($data_in[7]>$num_colloid)
	{
	    $num_colloid=$data_in[7];
	}
	#set all boundary particle to same sid/color.
	$data_in[7] = 1;	
	
    }
    elsif ($data_in[7] < 0 )
    {
	$data_in[7] = -1;
	$data_in[2] = 0;
	$data_in[3] = 0;
    }
    
    for($i=0;$i<$col;$i++)
    {
	$output[$n][$i]=$data_in[$i];
    }
    
    $n++;    
}

#change the sid/color of marked boundary particle

for($j=1;$j<$num_colloid;$j++)
{
    $output[$row[$j]][7]=2;
}

#output all the information
for($j=0;$j<$n;$j++)
{
    for($i=0;$i<$col;++$i)
    {
	#print $output[$j][$i], ' ',
	print file_out $output[$j][$i], ' '
    }
    #print "\n";
    print file_out "\n";
}


print $ARGV[1], "\n";

close(file_in);
close(file_out);
