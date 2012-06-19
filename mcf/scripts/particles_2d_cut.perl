#print $ARGV[0];
#print $ARGV[1];

$ARGV[1] = '>>' . $ARGV[1];
open file_in, $ARGV[0];
open file_out, $ARGV[1];

#$x_min=34;
#$x_max=46;
#$y_min=17;
#$y_max=23;
$x_min=14;
$x_max=26;
$y_min=14;
$y_max=26;

while($line = <file_in>)
{
 
    @data_in = split(' ', $line);

    if ($data_in[0]>  $x_min && 
	$data_in[0] < $x_max && 
	$data_in[1] > $y_min && 
	$data_in[1]<$y_max)
    {
	for($i=0;$i<=8;++$i)
	{
	    #print $data_in[$i], ' ',
	    print file_out $data_in[$i], ' '
	}
	#print "\n";
	print file_out "\n";
    }
}

print $ARGV[1], "\n";

close(file_in);
close(file_out);
