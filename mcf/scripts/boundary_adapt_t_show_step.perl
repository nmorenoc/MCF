
$dir=$ARGV[0];
open file_in, $dir."/mcf_boundary.dat";
$time=0.0;
$time=$ARGV[1];



while($line = <file_in>)
{
    
    @data = split(' ',$line);
    
    if ($data[1]>=$time)
    {
	print "step: ", $data[0], "\n";
	last;
    }
}
    
close(file_in);


