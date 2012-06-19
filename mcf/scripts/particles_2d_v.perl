#print $ARGV[0];
#print $ARGV[1];

$ARGV[1] = '>>' . $ARGV[1];
open file_in, $ARGV[0];
open file_out, $ARGV[1];

$scale = 6.27213e-6;

while($line = <file_in>)
{
 
    @data_in = split(' ', $line);
    
    
    if ( $data_in[7] == 0.0 ) 
    {
	
	$data_in[2]=sqrt($data_in[2]**2+$data_in[3]**2)/$scale;
    }
    else
    {
	$data_in[2]= -1;
    }
    
    
    for($i=0;$i<=2;++$i)
    {
	#print $data_in[$i], ' ',
	print file_out $data_in[$i], ' '
    }
    #print "\n";
    print file_out "\n";
    
}

print $ARGV[1], "\n";

close(file_in);
close(file_out);
