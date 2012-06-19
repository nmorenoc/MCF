open (fileinput1, "mcf_colloid01.dat");
open (fileinput2, "mcf_statistic.dat");
open (fileoutput, ">decay.dat");

while($line1 = <fileinput1>)
{
    $line2= <fileinput2>;

    @data1 = split(' ', $line1); 
    @data2 = split(' ', $line2); 
    
    $step = $data1[0];
    $t = $data1[1];
    $Ut = $data1[5];
    $Ucm = $data2[3];        
    
    print fileoutput  $step,' ',$t,' ', $Ut,' ',$Ucm, "\n";
}

close(fileinput1);
close(fileinput2);
close(fileoutput);

