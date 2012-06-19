open (fileinput, "particles_punto.dat");
open (fileoutput, ">particles_punto_sid.dat");

$n   = 0;
$col = 7;

while($line = <fileinput>)
{
    @data = split(' ', $line); 

    
    if($data[$col]<0) 
    {
	$data[$col] = 2;
	$n++;
    }   
    #print $t, ' ', $drag, ' ', "\n";
    print fileoutput $data[0],' ', $data[1],' ', $data[2],' ',$data[3],' ',$data[$col],"\n";
}

print "boundary particles    : ", $n, "\n";

close(fileinput);
close(fileoutput);

