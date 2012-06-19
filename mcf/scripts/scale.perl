#print $ARGV[0], "\n";

$filename_out = ">" . $ARGV[0] . ".scale";
#print $filename_out, "\n";

open file_in,  $ARGV[0];
open file_out, $filename_out;

$line = <file_in>;
@data_in=split(' ', $line);
#print $data_in[2];
$n_r    = $data_in[0];
$time_r = $data_in[2];
$speedup = $time_r/$data_in[2];
print file_out  $data_in[0], ' ', $data_in[1], ' ', $data_in[2], ' ', $speedup, "\n";

while($line=<file_in>)
{
    @data_in=split(' ', $line);
    #print $data_in[2];
    $speedup = $n_r*$time_r/$data_in[2];
    print file_out $data_in[0], ' ', $data_in[1], ' ', $data_in[2], ' ', $speedup, "\n"; 
}

close(file_in);
close(file_out);







