#print $ARGV[0];
#print $ARGV[1];

$ARGV[1] = '>>' . $ARGV[1];
open file_in, $ARGV[0];
open file_out, $ARGV[1];

$sid = 2;
$count = 0;

while($line = <file_in>)
{

    $count++;
    
    @data_in = split(' ', $line);
    
    if ($count == $sid ) 
    {
	for($i=0;$i<=12;++$i)
	{
	    #print $data_in[$i], ' ',
	    print file_out $data_in[$i], ' '
	}
	#print "\n";
	print file_out "\n";
	last;
    }
    
}

print $ARGV[1], "\n";

close(file_in);
close(file_out);
