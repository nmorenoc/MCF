#####################################################################
# Calculate local average velocity according to its distance
# to top and bottom boundaries.
#####################################################################


$dir_name=$ARGV[0];
$file_in_prefix="mcf_particles";
$file_out_prefix=$ARGV[1];

opendir(DIR, $dir_name) || die ("can not open directory");

print "processing directory : ", $dir_name, "\n";
@file_names_temp = readdir(DIR);
@file_names_in=sort @file_names_temp;
closedir(DIR);
print "number of files inside    : ", scalar(@file_names_in), "\n";

$num_file  = 0;
$step_start=50000;
$step_start=$ARGV[2];
$step_end  =9999999;
$step_end=$ARGV[3];

######################################
# Lx,Ly: box size.
# suppose x-direction periodic and
# y-direction wall.
######################################

$Lx=16.0;
$Ly=16.0;
$Lx=$ARGV[4];
$Ly=$ARGV[5];

#####################################################################
# interested gap in x-direction.
#####################################################################

$x_start= 0.0;
$x_end  = $Lx;
#print "x_start, x_end: ", $x_start, ' ',$x_end, "\n";


##################################################
# number of different resolutions.
##################################################
$num_res=1;

############
#resolution 0
############
$res[0]=80;
$res[0]=$ARGV[6];

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

########################
# Loop over files
#######################

$num_file=0;

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
	
	open (IN, $file_name);
	
	while ($line = <IN>)
	{
	    @data = split(' ', $line);
	    
	    $sx  = $data[0];
	    $sid = $data[7];
	    #print $sx, "\n";
	    
	    if( $sid==0 && $sx>=$x_start && $sx<=$x_end )
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
	
	$num_file++;
	close(IN);  
    }  
}

print "number of files processed: ", $num_file, "\n";


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



