open (fileinput, "conformation.dat");
open (fileoutput, ">verify.dat");

$dt  = 0.17361111E-02 ;

while($line = <fileinput>)
{
    @data = split(/ /, $line);
    
    $u1 = $data[7];
    $u2 = $data[8];
    $v1 = $data[9];
    $v2 = $data[10];
    
    $au1 = $data[11];
    $au2 = $data[12];
    $av1 = $data[13];
    $av2 = $data[14];
    
    $u1 = $u1 + 0.5 * dt *;

 print fileoutput $u1,' ', $u2,'(', $u1**2+$u2**2, ') * ', $v1, ' ', $v2, '(', $v1**2+$v2**2,') = ', $n, "\n"; 
}

close(fileinput);
close(fileoutput);
