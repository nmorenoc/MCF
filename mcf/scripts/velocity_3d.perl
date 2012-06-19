open (fileinput, "thin.dat");
open (fileoutput,">v.dat");

$n = 0;
$v1 = 0;
$v2 = 0;
$v3 = 0;


while($line = <fileinput>)
{
 @data = split(/ /, $line);
 $v1 = $v1 + $data[3];
 $v2 = $v2 + $data[4];
 $v3 = $v3 + $data[5];
 $n = $n +1;

 if (@data[0] == '') 
  {
  	$v1 = $v1 /$n;
	$v2 = $v2 /$n;
	$v3 = $v3 /$n;
#	print $v1,' ', $v2,' ',$v3, "\n";	
	print fileoutput  $v1,' ',$v2,' ',$v3, "\n";	
	$v1 = 0;
	$v2 = 0;
	$v3 = 0;
	$n = 0;	
   }

}

close(fileinput);
close(fileoutput);
