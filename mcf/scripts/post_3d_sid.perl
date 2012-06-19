#print $ARGV[0];
#print $ARGV[1];

$ARGV[1] = '>>' . $ARGV[1];
open file_in, $ARGV[0];
open file_out, $ARGV[1];

print $ARGV[1];

while($line = <file_in>)
{
 
    @data_in = split(' ', $line);
    
    if ($data_in[9] > 0)
    {
	#$data_in[9] = 1.0;
    }
    elsif ($data_in[9] < 0 ) 
    {
	$data_in[9] = -1.0;
	#next;
    }
    elsif($data_in[9] == 0 )
    {
	next;
    }
	
    
    for($i=0;$i<=12;++$i)
    {
	#print $data_in[$i], ' ',
	print file_out $data_in[$i], ' '
    }
    #print "\n";
    print file_out "\n";
}

close(file_in);
close(file_out);
