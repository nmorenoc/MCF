#print $ARGV[0];
#print $ARGV[1];

$ARGV[1] = '>>' . $ARGV[1];
open file_in, $ARGV[0];
open file_out, $ARGV[1];


$n=0.0;

while($line = <file_in>)
{
 
    @data_in = split(' ', $line);

    for($i=0;$i<=5;++$i)
    {
	#print $data_in[$i], ' ',
	print file_out $data_in[$i], ' ';
	
    }
    $n=$n+0.1;
    print file_out $n;
    #print "\n";
    print file_out "\n";
}

print $ARGV[1], "\n";

close(file_in);
close(file_out);
