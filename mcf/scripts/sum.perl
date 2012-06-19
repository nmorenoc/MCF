open (fileinput, "sum.dat");

$num_tot = 0;
$n = 0;

while($line = <fileinput>)
{
    @data = split(' ', $line); 
    
    $num = $data[6];
    print $data[6], "\n";
    $n = $n + 1;
    $num_tot = $num_tot + $num;
    
}

print "n       : ", $n, "\n";
print "num_tot : ", $num_tot, "\n";

close(fileinput);

