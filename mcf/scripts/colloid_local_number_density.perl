#############################################################
# Calculate local volume fraction of colloid according to 
# its distance to top and bottom boundaries.
#############################################################

$dir_name= $ARGV[0];
$file_in_prefix="mcf_colloid";
$file_out_prefix= $ARGV[1];

$R=1.0;
$pi=3.1415926;

opendir(DIR, $dir_name) || die ("can not open directory");
print "processing directory : ", $dir_name, "\n";
@file_names_temp = readdir(DIR);
@file_names_in=sort @file_names_temp;
closedir(DIR);


$num_file  = 0;
$step_start=300000;
$step_start=$ARGV[2];
$step_end  =99999999;
$step_end  =$ARGV[3];

$Lx=16.0;
$Ly=16.0;
$Lx=$ARGV[4];
$Ly=$ARGV[5];

#############################################################
# number of different resolutions.
#############################################################
$num_res=3;

#############################################################
#the 0th resolution 
#############################################################

$res[0]=40;
$res[0]=$ARGV[6];

$h[0]=$Ly/$res[0];

print "number of resolution : ", $num_res, "\n";
print "reslutions : ", $res[0], " ";

#############################################################
# further resolution is two times than previous one.
# Note it is questionable to use higher resolution
# than SPH/SDPD resolution.
#############################################################
for ($j=1;$j<$num_res;$j++)
{
    $res[$j]=$res[$j-1]*2;
    $h[$j]=$Ly/$res[$j];
    print $res[$j], " ";
}
print "\n";


#############################################################
# initialize mid location y and number counter of each layer.
#############################################################

for ($j=0;$j<$num_res;$j++)
{
    for ($k=0; $k<$res[$j]; $k++)
    {
	$y[$k][$j] = $k*$h[$j]+ $h[$j]/2.0;
	$num[$k][$j] = 0.0;
    }
}

#############################################################
# Loop over all files according to starting and ending files. 
#############################################################

$num_file=0;

print "starting step : ", $step_start;
print "ending  step  : ", $step_end;
 
$f_start = $file_in_prefix . stepstring($step_start). ".out";
$f_end = $file_in_prefix . stepstring($step_end). ".out";

print "starting file : ", $f_start, "\n";
print "ending file   : ", $f_end, "\n";

foreach $f (@file_names_in)
{
    if (( $f ge $f_start) && ( $f lt $f_end ) )
    {
	$file_name = $dir_name . $f;
	print "processing file : ", $file_name, "\n";
	
#############################################################
# reset to zero 
#############################################################
	
	open (IN, $file_name);
	
	$num_particle=0;

	while ($line = <IN>)
	{
	    @data = split(' ', $line);
	    
	    $y = $data[1];
	    
	    for ($j=0; $j<$num_res; $j++)
	    {
		#print "y/h : ", $y/$h, "\n";
		$num[$y/$h[$j]][$j] += 1.0;
	    }
	    
	    $num_particle++;
	    
	}
	$num_file ++;
	
	close(IN);
	
    }
    
}

print "number of files processed: ", $num_file, "\n";

if ($num_file > 0 ) 
{
    for ($j=0; $j<$num_res; $j++)
    {
	for ($k=0;$k<=$res[$j];$k++)
	{
	    $num[$k][$j]=$num[$k][$j]*$res[$j]/$num_file/$num_particle;
	}
    }
}

for ($j=0; $j<$num_res; $j++)
{
    $file_out_name = ">".$file_out_prefix."_".$h[$j].".dat";
    open(OUT, $file_out_name);

    for ($k=0;$k<$res[$j];$k++)
    {
	print OUT $y[$k][$j],' ', $num[$k][$j], "\n";
    }
    print "writing file : ", $file_out_name, "\n";
    close(OUT);    
}


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

