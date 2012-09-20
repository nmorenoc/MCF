
$filein_name=$ARGV[0];

open (fileinput, $filein_name);


$num_line = 0;
$num_dt   = 0;
$freq     = 1;

$stress_av=0;
$rdf_av   =0;

$stress2_av= 0;
$rdf2_av=0;

while($line = <fileinput>)
{
    $num_line = $num_line + 1;

    if($num_line%$freq==0) 
    {
	$num_dt = $num_dt + 1;
	
	@line_data = split(' ', $line); 
	push @data,[$line_data[0],$line_data[1],$line_data[2],$line_data[3]];

	$stress_av += $line_data[2];
	$rdf_av += $line_data[3];

	$stress2_av += $line_data[2]**2;
	$rdf2_av += $line_data[3]**2;
	
    }
}

#average
$stress_av /= $num_dt;
$rdf_av /= $num_dt;
$stress2_av /= $num_dt;
$rdf2_av /= $num_dt;

#standard variance
$stress_var=sqrt($stress2_av-$stress_av**2);
$rdf_var=sqrt($rdf2_av-$rdf_av**2);


for( $n_dt = 0; $n_dt<$num_dt; $n_dt++)
{
    $v_prod = 0.0;
    $count  = 0.0;
    
    for ($s0=0; $s0<$num_dt-$n_dt; $s0++)
    {
	$s1=$s0+$n_dt;
	$dot_prod = (($data[$s0]->[2])-$stress_av)*(($data[$s1]->[3])-$rdf_av);
	$v_prod = $v_prod + $dot_prod;
	$count = $count + 1.0 ;	
    }

    $v_prod = $v_prod/$count;
    print  $n_dt, ' ', $v_prod/$stress_var/$rdf_var, "\n";
}

close(fileinput);

