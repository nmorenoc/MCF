open file_in1, "x.out";
open file_in2, "f.out";
open file_in3, "id.out";
open file_out, ">data.dat";


$n = 0;

while( $line_in1 = <file_in1> )
{
    
    $line_in2 = <file_in2> ;
    $line_in3 = <file_in3> ;

#    print $line_in1, $line_in2, $line_in3, "\n";
    
    @data_in1 = split(' ',$line_in1);
    @data_in2 = split(' ',$line_in2);
    @data_in3 = split(' ',$line_in3);

#   print $data_in[0], "\n";
    if ( $data_in3[1] == 0) 
    {
	$data_in2[0] = 0;
	$data_in2[1] = 0;
	    
    }
    print file_out $data_in1[0], ' ', $data_in1[1], ' ', $data_in2[0], ' ', $data_in2[1], ' ', $data_in3[0], ' ', $data_in3[1],"\n"; 
    $n=$n+1;
}
    
close(file_in1);
close(file_in2);
close(file_in3);
(file_out);

