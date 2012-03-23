#############################################################
# Calculate local volume fraction of colloid according to 
# its distance to top and bottom boundaries.
#############################################################

#$dir_name="../C0589_32_16_repuln1_2C/";
#$file_in_prefix="mcf_colloid";
#$file_out_prefix="lnd_C0589_32_16";

$dir_name_prefix= $ARGV[0];
$dir_num=$ARGV[1];
$file_in_prefix="mcf_colloid";
$file_out_prefix= $ARGV[2];
$Lx=16.0;
$Ly=16.0;
$R=1.0;
$pi=3.1415926;

$num_file  = 0;

$step_start=300000;
$step_start=$ARGV[3];
$step_end  =99999999;
$step_end  =$ARGV[4];

#############################################################
# number of different resolutions.
#############################################################
$num_res=6;

#############################################################
#the 0th resolution 
#############################################################

$res[0]=40;
#$res_step=40;
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
# loop each folder given.
#############################################################


for ($fd=0;$fd<=$dir_num;$fd++)
{
    $dir_name=$dir_name_prefix."_".$fd."/";
    print "processing directory : ", $dir_name, "\n";
    opendir(DIR, $dir_name) || die ("can not open directory");
    @file_names_temp = readdir(DIR);
    @file_names_in=sort @file_names_temp;
    closedir(DIR);
    
    
#############################################################
# Loop over all files according to starting and ending files. 
#############################################################
    
    $num_file[$fd]=0;
    
    $f_start = $file_in_prefix . stepstring($step_start). ".out";
    $f_end = $file_in_prefix . stepstring($step_end). ".out";
    
    print "starting file : ", $f_start, "\n";
    print "ending file : ", $f_end, "\n";
    
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
	    
	    $num_file[$fd] ++;
	    
	    close(IN);
	    
	}
	
    } # each file
}#each folder


$num_file_tot=0;
for ($fd=0;$fd<=$dir_num;$fd++)
{
    $num_file_tot+=$num_file[$fd];
    print "number of files processed: ", $num_file[$fd], "\n";    
}
print "total number of files processed: ", $num_file_tot, "\n";

if ($num_file_tot > 0 ) 
{
    for ($j=0; $j<$num_res; $j++)
    {
	for ($k=0;$k<=$res[$j];$k++)
	{
	    $num[$k][$j]=$num[$k][$j]*$res[$j]/$num_file_tot/$num_particle;
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

