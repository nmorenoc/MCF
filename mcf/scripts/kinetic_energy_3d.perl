open (fileinput, "mcf_particles00150000.out");


$k1 = 0.0;
$k2 = 0.0;
$k3 = 0.0;


while($line = <fileinput>)
{
 @data = split(' ', $line);
 
 $v1 = $data[3];
 $v2 = $data[4];
 $v3 = $data[5];
 $m  = $data[7];

 #print @data, "\n";
 #print ' ', $m, ' ', $v1, ' ', $v2, ' ', $v3, "\n";
 #last;

 $k1 = $k1 + 0.5*$m*$v1**2;
 $k2 = $k2 + 0.5*$m*$v2**2;
 $k3 = $k3 + 0.5*$m*$v3**2;

}

print 'k1,k2,k3 : ', $k1,' ', $k2, ' ', $k3,  "\n";
print 'k_tot    : ' ,$k1+$k2+$k3, "\n";

close(fileinput);
