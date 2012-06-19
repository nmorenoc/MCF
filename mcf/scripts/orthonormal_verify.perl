open (fileinput, "conformation.dat");
open (fileoutput, ">verify.dat");

while($line = <fileinput>)
{
 @data = split(/ /, $line);
 $u1 = $data[3];
 $u2 = $data[4];
 $v1 = $data[5];
 $v2 = $data[6];

 $n = $u1*$v1 + $u2*$v2;

 print fileoutput $u1,' ', $u2,'(', $u1**2+$u2**2, ') * ', $v1, ' ', $v2, '(', $v1**2+$v2**2,') = ', $n, "\n"; 
}

close(fileinput);
close(fileoutput);
