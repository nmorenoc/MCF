#############################################################
# Count averaged number of neigbors at different
# distances in a suspension.
#############################################################
use Math::Trig;


$dir_name=$ARGV[0];
$file_in_prefix="mcf_colloid";
$file_out_prefix=$ARGV[1];

opendir(DIR, $dir_name) || 
die   ("can not open directory   :".$dir_name.", as it does not exist !\n");
print "processing directory      : ", $dir_name, "\n";

@file_names_temp = readdir(DIR);
@file_names_in=sort @file_names_temp;
closedir(DIR);

$num_file_tot  = scalar(@file_names_in);
print "number of files inside    : ", $num_file_tot, "\n";


$pi=3.1415926;

$step_start=80000;
$step_start=$ARGV[2];
$step_end  =99999999;
$step_end=$ARGV[3];

######################################
# Lx,Ly: box size.
# V    : volume 3D / area 2D.
# gap  : gap considered away from the walls.
# suppose x-direction periodic and
# y-direction wall.
# dt   : time step
# R    : radius
######################################

$Lx=$ARGV[4];
$Ly=$ARGV[5];
$gap=$ARGV[6];
$dt=1.329787234042553e-3;
$dt=$ARGV[7];
$R=1.0;

################################################################
#num_p: number of colloid in the file.
#num_o: number of colloid considered as origin, 
#       for which we look for neighbors.
#num_t: number of all including the ones near wall and periodic
################################################################
$num_p=0;
$num_o=0;
$num_t=0;

#############################################################
# minimum/maximum distrance.
#############################################################
$r_min=2.0;
$r_max=2.5;

#############################################################
# number of different resolutions in distance.
#############################################################
$num_dr=11;

print "minimum and maximum distance, number of layers: ", 
    $r_min, ' ', $r_max, ' ', $num_dr,;
print "\n";

#############################################################
# distance of each layer.
#############################################################

print "distance at each layer: ";

for ($i=0;$i<$num_dr;$i++)
{
    $r[$i]=$r_min+$i*($r_max-$r_min)/($num_dr-1);
    print $r[$i], ' ';
}
print "\n";



#############################################################
# initialize number counter of each portion.
# first dimension is  for step, time and different distances.
# second dimension is for each time step.
#############################################################

for ($i=0; $i<$num_file_tot; $i++)
{
    for ($j=0; $j<$num_dr+2; $j++)
    {
	$num[$j][$i] = 0.0;	
    }
}

#############################################################
# Loop over all files according to starting and ending files.
#############################################################

$num_file=0;

$f_start = $file_in_prefix . stepstring($step_start). ".out";
$f_end = $file_in_prefix . stepstring($step_end). ".out";
print "starting file   : ", $f_start, "\n";
print "ending file     : ", $f_end, "\n";


foreach $f (@file_names_in)
{
    #print "checking file : ", $f, "\n";
    
    if (( $f ge $f_start) && ( $f le $f_end ) )
    {
	$step=stringstep($f);
	$file_name = $dir_name . $f;
	print "processing file : ", $file_name, "\n";
	
#############################################################
# Open a colloid file.
# reset counter to zero. 
# read each colloid position.
#############################################################
	
	open (IN, $file_name);	
	
	$num_p=0;
	while ($line = <IN>)
	{
	    @data = split(' ', $line);
	    
	    $x[$num_p] = $data[0];
	    $y[$num_p] = $data[1];
	    $num_p++;
	}
	close(IN);
	#print "num_p : ", $num_p, "\n";
	
#############################################################
# reset counter to zero. 
# record the origin ones, which are inside simulation box and
# away from walls.
# record all ones, assuming x direction is periodic boundary.
#############################################################
	$num_o=0;
	$num_t=0;
	
	for ($p=0;$p<$num_p;$p++)
	{
	    if ( ( $y[$p]>$gap ) && ( $y[$p]< $Ly-$gap) )
	    {
		$x_o[$num_o] = $x[$p];
		$y_o[$num_o] = $y[$p];
		$num_o++;
	    }
	    $x_t[$num_t]=$x[$p];
	    $y_t[$num_t]=$y[$p];
	    $num_t++;
	    $x_t[$num_t]=$x[$p]-$Lx;
	    $y_t[$num_t]=$y[$p];
	    $num_t++;
	    $x_t[$num_t]=$x[$p]+$Lx;
	    $y_t[$num_t]=$y[$p];
	    $num_t++;
	}
        #print "num_o : ", $num_o, "\n";
	#print "num_t : ", $num_t, "\n";

#############################################################
# calculate number of neigbours according to different
# distances.
# periodic images have to be considered.
#############################################################

	#####################################
	# record step, time first
	#####################################
	
	$num[0][$num_file] = $step;
	$num[1][$num_file] = $step*$dt;

	for($p=0;$p<$num_o; $p++)
	{  
	    for($q=0;$q<$num_t;$q++)
	    {
		$x_12 = $x_o[$p]-$x_t[$q];
		$y_12 = $y_o[$p]-$y_t[$q];
		$r_12 = sqrt($x_12**2+$y_12**2);
		
		#print "r_12: ", $r_12, "\n";
		#####################################
                # Count neighbors, actually itself
                # is also included.
		#####################################

		for ($i=0;$i<$num_dr;$i++)
		{
		    if( $r_12 <= $r[$i] )
		    {
			$num[2+$i][$num_file] ++;
		    }
		    
		} # i < num_dr
	    } # q < num_t
	} # p < num_o
	
#############################################################
# Average over num_o particles and exclude itself.
#############################################################
	
	for ($i=0;$i<$num_dr;$i++)
	{	    
	    $num[2+$i][$num_file] /= $num_o; 
	    $num[2+$i][$num_file] -= 1;
	} # i < num_dr
	
#############################################################	
# increase counter of files
#############################################################
	$num_file ++;
	
    } # if file is in consideration.
    
}# folder

print "number of files processed: ", $num_file, "\n";



$file_out_name = ">".$file_out_prefix."_".$gap.".dat";
open(OUT, $file_out_name);

for ($i=0;$i<$num_file;$i++)
{
    for ($j=0; $j<$num_dr+2; $j++)
    {
	print OUT $num[$j][$i], ' ';
    }
    
    print OUT "\n";    
}
    
print "writing output file : ", $file_out_name, "\n";
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

sub stringstep
{

    $string=$_[0];
    
    $sub_string = substr $string, 11, 8;
    
    #print $sub_string, "\n";
    $step = int($sub_string);
    #print $step, "\n";
    return $step;
}

