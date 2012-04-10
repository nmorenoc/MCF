$dir_name= $ARGV[0];
$file_in_prefix="mcf_colloid";
open(OUT, ">mcf_colloid0001.dat");

$step_start=0;
$step_start=$ARGV[1];
$step_end  =99999999;
$step_end  =$ARGV[2];

$dt=1.050696086157079E-003;
$freq=10;
$num_col=18;

opendir(DIR, $dir_name) || die ("can not open directory");
print "processing directory : ", $dir_name, "\n";
@file_names_temp = readdir(DIR);
@file_names_in=sort @file_names_temp;
closedir(DIR);

$ss=0;

$f_start = $file_in_prefix . stepstring($step_start). ".out";
$f_end = $file_in_prefix . stepstring($step_end). ".out";

print "f_start: ", $f_start, "\n";
print "f_end:   ", $f_end, "\n";
foreach $f (@file_names_in)
{
    if (( $f ge $f_start) && ( $f lt $f_end ) )
    {
	#$file_name = $dir_name . $f;
	$file_name = $f;
	print "processing file : ", $file_name, "\n";

	open(IN,$file_name);
	
	while ($line = <IN>)
	{
	    @data = split(' ', $line);
	}
	
	#print @data[0..10], "\n";

	$time=$ss*$dt;
	
	
	print OUT $ss, ' ', $time, ' ';
	for ($i=0;$i<$num_col; $i++)
	{
	    print OUT $data[$i], ' ';
	}
	print OUT "\n";
	
	$ss+=$freq;

    }
}

close(OUT);

sub stepstring
{

    $step=$_[0];
    $length_step=length($step);

    $string = $step;
    for ($i=$length_step;$i<8;$i++)
    {
	$string = "0" . $string;
    }

    return $string;
}
