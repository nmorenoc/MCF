open file_in, "msd.dat";
open file_t, ">msd_colloid01.dat";


$dt = 0.125E-04;
$n  = 0;

while($line_in = <file_in>)
{
    
 #   print $line_in, "\n";
    
    @data_in = split(' ',$line_in);
    
#   print $data_in[0], "\n";
    print file_t  $dt*$n, ' ',$data_in[0], "\n"; 
    $n=$n+1;
}
    
close(file_t);
close(file_in);

