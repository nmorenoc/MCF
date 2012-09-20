#####################################################################
# Calculate local average velocity according to its distance
# to top and bottom boundaries.
#####################################################################


$dir_prefix=$ARGV[0];
$num_dir=$ARGV[1];
$file_in_prefix="mcf_colloid";
$file_out_prefix=$ARGV[2];

$step_start=50000;
$step_start=$ARGV[3];
$step_end  =9999999;
$step_end=$ARGV[4];

######################################
# Lx,Ly: box size.
# suppose x-direction periodic and
# y-direction wall.
######################################

$Lx=16.0;
$Ly=16.0;
$Lx=$ARGV[5];
$Ly=$ARGV[6];

#####################################################################
# interested gap in x-direction.
#####################################################################

$x_start= 0.0;
$x_end  = $Lx;
#print "x_start, x_end: ", $x_start, ' ',$x_end, "\n";


##################################################
# number of different resolutions.
##################################################
$num_res=3;

############
#resolution 0
############
$res[0]=40;
$h[0]=$Ly/$res[0];

print "number of resolution : ", $num_res, "\n";
print "res : ", $h[0], ' ',$res[0], "; ";

############
#other resolutions
############
for ($j=1;$j<$num_res;$j++)
{
    $res[$j]=$res[$j-1]*2.0;
    $h[$j]=$Ly/$res[$j];
    print $h[$j], ' ', $res[$j], "; ";
}
print "\n";

########################
# initialize y and v
########################


for ($j=0;$j<$num_res;$j++)
{
    for ($k=0; $k<$res[$j]; $k++)
    {
	$y[$k][$j] = $k*$h[$j]+ $h[$j]/2.0;
	$v[$k][$j] = 0.0;
	$n[$k][$j] = 0.0
    }
}



print "starting step : ", $step_start;
print "ending  step  : ", $step_end;
  
$f_start = $file_in_prefix . stepstring($step_start). ".out";
$f_end = $file_in_prefix . stepstring($step_end). ".out";
print "starting file : ", $f_start, "\n";
print "ending file : ", $f_end, "\n";

###################################################
#loop each folder given.
###################################################

for ($fd=1;$fd<=$num_dir;$fd++)
{
    $dir_name=$dir_prefix."_".$fd."/";
    opendir(DIR, $dir_name) || die ("can not open directory");
    print "processing directory : ", $dir_name, "\n";

    @file_names_temp = readdir(DIR);
    @file_names_in=sort @file_names_temp;
    print "number of files inside    : ", scalar(@file_names_in), "\n";
    closedir(DIR);
    
    
########################
# Loop over files
#######################
    
    $num_file[$fd]=0;

    foreach $f (@file_names_in)
    {
	if (( $f ge $f_start) && ( $f lt $f_end ) )
	{
	    $file_name = $dir_name . $f;
	    print "processing file : ", $file_name, "\n";
	    
	    open (IN, $file_name);
	    
	    while ($line = <IN>)
	    {
		@data = split(' ', $line);
		
		$sx  = $data[0];
		#print $sx, "\n";
		
		if( $sx>=$x_start && $sx<=$x_end )
		{
		    $sy = $data[1];
		    $vx = $data[2];
		    
		    for ($j=0; $j<$num_res; $j++)
		    {
			#print "sy/$h : ", $sy/$h[$j], "\n";
			$v[$sy/$h[$j]][$j] += $vx;
			$n[$sy/$h[$j]][$j] ++;
		    }
		}
	    }
	    
	    $num_file[$fd]++;
	    close(IN);  
	}  
    }# each file
}#each folder

$num_file_tot=0;
for ($fd=1;$fd<=$num_dir;$fd++)
{
    $num_file_tot+=$num_file[$fd];
    print "number of files processed: ", $num_file[$fd], "\n";
}

print "total number of files processed: ", $num_file_tot, "\n";

for ($j=0; $j<$num_res; $j++)
{
    
    $file_out_name = ">".$file_out_prefix."_".$h[$j].".dat";
    open(OUT, $file_out_name);
    
    for ($k=0; $k<$res[$j]; $k++)
    {
	if ( $n[$k][$j] == 0 )
	{
	    $v[$k][$j] = 0.0;
	}
	else
	{
	    $v[$k][$j] = $v[$k][$j]/$n[$k][$j];
	}
	print OUT $y[$k][$j],' ', $v[$k][$j], "\n";
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



