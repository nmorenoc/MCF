open file_in, "mcf_boundary.dat";

$time=0.0;

$time=$ARGV[0];



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


