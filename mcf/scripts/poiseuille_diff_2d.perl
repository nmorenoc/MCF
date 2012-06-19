open file_n, "thin.dat";
open file_a, "poiseuille_b.dat";

$n = 0;
$L1 = 0.0;
$L1_t = 0.0;
$L2 = 0.0;
$abs_diff_aver = 0.0;
$abs_diff_max = 0.0;
$abs_diff_min = 1000.0;

$diff_aver = 0.0;
$diff_max = 0.0;
$diff_min = 1000.0;

while($line_n = <file_n> and $line_a = <file_a>)
{

#    print $line_n, "\n";
#    print $line_a, "\n";

    @data_n = split(' ', $line_n);
    @data_a = split(' ', $line_a);
    
    #print $data_n[0], "\n";
    #print $data_a[1], "\n";
    
    $L1 =  $L1+abs($data_n[2] - $data_a[1]);
    $L1_t = $L1_t+$data_a[1];
    $L2 =  $L2 + ($data_n[2] - $data_a[1])**2;
    $abs_diff_aver = $abs_diff_aver + ($data_n[2] - $data_a[1]) / $data_a[1];
    if(  abs(($data_n[2] - $data_a[1]))/$data_a[1] > abs($abs_diff_max))
    {
	$abs_diff_max = abs(($data_n[2] - $data_a[1])/$data_a[1]);
    }    
    if(  abs(($data_n[2] - $data_a[1]))/$data_a[1] < abs($abs_diff_min))
    {
	$abs_diff_min = abs(($data_n[2] - $data_a[1])/$data_a[1]);
    }  

    $diff_aver = $diff_aver + ($data_n[2] - $data_a[1]) / $data_a[1];
    if(  ($data_n[2] - $data_a[1])/$data_a[1] > $diff_max)
    {
	$diff_max = ($data_n[2] - $data_a[1])/$data_a[1];
    }    
    if(  abs(($data_n[2] - $data_a[1]))/$data_a[1] < $diff_min)
    {
	$diff_min = ($data_n[2] - $data_a[1])/$data_a[1];
    }  


  
    $n = $n + 1;	
}
$L2 = sqrt($L2);
$abs_diff_aver = $abs_diff_aver / $n;
$diff_aver = $diff_aver / $n;

print "n:", $n, "\n";
print "L1:", $L1, "\n";
print "L1-norm:", $L1/$L1_t,"\n";
print "L2:", $L2, "\n";
print "abs_diff_aver:", $abs_diff_aver, "\n";
print "abs_diff_max:", $abs_diff_max, "\n";
print "abs_diff_min:", $abs_diff_min, "\n";

print "diff_aver:", $diff_aver, "\n";
print "diff_max:", $diff_max, "\n";
print "diff_min:", $diff_min, "\n";

close(file_n);
close(file_a);

