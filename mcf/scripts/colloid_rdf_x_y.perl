#############################################################
# Calculate radial distribution function of colloid pair
# according to its distance and relative angle and
# map on a 2D grid;
#############################################################


$dir_name=$ARGV[0];
$file_in_prefix="mcf_colloid";
$file_out_prefix=$ARGV[1];

opendir(DIR, $dir_name) || die ("can not open directory");
print "processing directory : ", $dir_name, "\n";

@file_names_temp = readdir(DIR);
@file_names_in=sort @file_names_temp;
closedir(DIR);
print "number of files inside     : ", scalar(@file_names_in), "\n";


$num_file  = 0;
$step_start=8000;
$step_start=$ARGV[2];
$step_end  =99999999;
$step_end=$ARGV[3];

#############################################################
# Box size of simulation
# gap: remove colloids which have distance smaller than gap.
#############################################################
$Lx=$ARGV[4];
$Ly=$ARGV[5];
$V=$Lx*$Ly;
$gap=$ARGV[6];
$R=1.0;

#############################################################
# minimum/maximum distrance
# minimum angle=0; maximum angle=2*pi.
#############################################################
$r_min=-5.0;
$r_max=5.0;

#############################################################
# number of different resolutions in distance.
#############################################################
$num_dr=3;

#############################################################
#the 0th resolution, i.e., number of points.
#############################################################

$r_num[0]=50;
$r_dr[0]=($r_max-$r_min)/$r_num[0];

print "number of resolution : ", $num_dr, "\n";

#############################################################
#other resolutions
#############################################################

for ($j=1;$j<$num_dr;$j++)
{
    $r_num[$j]=$r_num[$j-1]*2;
    $r_dr[$j] =($r_max-$r_min)/$r_num[$j];
}


print "reslutions (x_num,y_num,dx,dy;): ";

for ($i=0;$i<$num_dr;$i++)
{
    print  $r_num[$i], ' ',$r_num[$i], ' ', $r_dr[$i], ' ',$r_dr[$i],"; ";    
}

print "\n";


#############################################################
# initialize mid location r, t.
# as x and y dimension have the same resolution,
# use one dimension array to record its mid location,
# use second dimension for different resolution.
#############################################################

for ($j=0;$j<$num_dr;$j++)
{
    for ($k=0; $k<$r_num[$j]; $k++)
    {
	$r[$k][$j] = $r_min+$k*$r_dr[$j]+ $r_dr[$j]/2.0;
	#print $r[$k][$j], "\n";
    }
}

#############################################################
# initialize number counter of each portion.
#############################################################

for ($i=0; $i<$num_dr;$i++)
{
    for ($j=0; $j<$r_num[$i];$j++)
    {
	for ($k=0; $k<$r_num[$i];$k++)
	{
	    $num[$j][$k][$i] = 0.0;
	    
	}
    }
}

#############################################################
# Loop over all files according to starting and ending files.
#############################################################

$num_file=0;

$f_start = $file_in_prefix . stepstring($step_start). ".out";
$f_end = $file_in_prefix . stepstring($step_end). ".out";
print "starting file : ", $f_start, "\n";
print "ending file   : ", $f_end, "\n";

################################################################
#num_p: number of colloid in the file.
#num_o: number of colloid considered as origin, 
#       for which we look for neighbors.
#num_t: number of all including the ones near wall and periodic
################################################################
$num_p=0;
$num_o=0;
$num_t=0;

foreach $f (@file_names_in)
{
    #print "checking file : ", $f, "\n";
    
    if (( $f ge $f_start) && ( $f le $f_end ) )
    {
	$file_name = $dir_name . $f;
	print "processing file : ", $file_name, "\n";
	
#############################################################
# Open a colloid file.
# reset counter to zero. 
# ### the last simulation file may be incomplete
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
	
	###################################
        # get number denisty.
        ###################################
	$num_density=$num_p/$V;
	$num_density_r=$V/$num_p;
	#print "num density : ", $num_density, "\n";

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
# calculate pair-wise correlation.
# i.e., distance and relative angle.
# periodic images have to be considered.
#############################################################
	
	$coeff = $num_density_r/$num_o/2.0;
	
	for($p=0;$p<$num_o; $p++)
	{
	    for($q=$0;$q<$num_t;$q++)
	    {
		
		$x_12 = $x_o[$p]-$x_t[$q];
		$y_12 = $y_o[$p]-$y_t[$q];
		
		#print "x_12,y_12: ", $x_12, ' ', $y_12, "\n";
		
		for ($i=0;$i<$num_dr;$i++)
		{
		    ###########################################
                    # Get relative position.
		    ###########################################
		    
		    $x_idx=($x_12-$r_min)/$r_dr[$i];
		    $y_idx=($y_12-$r_min)/$r_dr[$i];

		    ###########################################
                    # exclude two particles too far.
                    ###########################################
		    
		    if( $x_idx<0 || $x_idx >= $r_num[$i] )
		    {
			next;
		    }
		    if( $y_idx<0 || $y_idx >= $r_num[$i] )
		    {
			next;
		    }
		   
		    #print "x_idx, y_idx : ", $x_idx, ' ' , $y_idx, "\n";
		    $num[$x_idx][$y_idx][$i] += ($coeff/$r_dr[$i]/$r_dr[$i]);
		    
		    ###########################################
                    # Get reverse relative position.
		    ###########################################

		    $x_idx=(-$x_12-$r_min)/$r_dr[$i];
		    $y_idx=(-$y_12-$r_min)/$r_dr[$i];
		    
		    ###########################################
                    # exclude two particles too far.
                    ###########################################
		    
		    if( $x_idx<0 || $x_idx >= $r_num[$i] )
		    {
			next;
		    }
		    if( $y_idx<0 || $y_idx >= $r_num[$i] )
		    {
			next;
		    }
		    
		    #print "x_idx, y_idx : ", $x_idx, ' ' , $y_idx, "\n";
		    ############################################
                    # Add one contribution.
		    ############################################
		    $num[$x_idx][$y_idx][$i] += ($coeff/$r_dr[$i]/$r_dr[$i]);
		    
		} # i < num_dr
	    } # q < num_t
	} # p < num_o
	
#############################################################	
# increase counter of files
#############################################################
	$num_file ++;
	
    } # if file is in consideration.
    
}# folder

#print "number of particles: ", $num_p,"\n";
print "number of files processed: ", $num_file, "\n";

if ( $num_file > 0) 
{
    for ($i=0; $i<$num_dr; $i++)
    {	
	for ($j=0;$j<$r_num[$i];$j++)
	{
	    for ($k=0;$k<$r_num[$i];$k++)
	    {
		$num[$j][$k][$i] /= $num_file;
	    }
	}
	
    }
}

for ($i=0; $i<$num_dr; $i++)
{
    $file_out_name = ">".$file_out_prefix."_".$gap."_".$r_dr[$i].".dat";
    open(OUT, $file_out_name);
    
    for ($j=0;$j<$r_num[$i];$j++)
    {
	for ($k=0;$k<$r_num[$i];$k++)
	{
	    $dij2 = $r[$j][$i]**2 + $r[$k][$i]**2;
	    
	    if ($dij2<=$R**2)
	    {
		$num[$j][$k][$i] = -1;
	    }
	    print OUT $r[$j][$i],' ', $r[$k][$i], ' ', $num[$j][$k][$i], "\n";
	}
	print OUT "\n";
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

