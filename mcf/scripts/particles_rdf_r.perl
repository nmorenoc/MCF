#############################################################
# Calculate radial distribution function of particles pair
# according to its radial distance, it is averaged over
# orientation angle theta.
#############################################################
use Math::Trig;


$dir_name=$ARGV[0];
$file_in_prefix="mcf_particles";
$file_out_prefix=$ARGV[1];

opendir(DIR, $dir_name) || 
die ("can not open directory, it does not exist !\n");

print "processing directory      : ", $dir_name, "\n";
@file_names_temp = readdir(DIR);
@file_names_in=sort @file_names_temp;
closedir(DIR);
print "number of files inside    : ", scalar(@file_names_in), "\n";


$pi=3.1415926;
$num_file  = 0;
$step_start=22560;
$step_start=$ARGV[2];
#$step_end  =200000;
$step_end  =99999999;

######################################
# Lx,Ly: box size.
# V    : volume 3D / area 2D.
# gap  : gap considered away from the walls.
# suppose x-direction periodic and
# y-direction wall.
######################################

$Lx=4.0;
$y_min_t=-1.0;
$y_max_t=5.0;
$y_min=0.0;
$y_max=4.0;
$Ly=$y_max-$y_min;
$V=$Lx*$Ly;
$gap=0.2;
$R=0.2;
$V_eff=$Lx*($Ly-2.0*($gap-$R));


################################################################
#num_p:  number of colloid in the file.
#num_o: number of colloid considered as origin, 
#       for which we look for neighbors.
#num_t: number of all including the ones near wall and periodic
################################################################
$num_p=0;
$num_o=0;
$num_t=0;

#############################################################
# minimum/maximum distrance
#############################################################
$r_min=0.2;
$r_max=2.0;

#############################################################
# number of different resolutions in distance.
#############################################################
$num_dr=4;

#############################################################
#the 0th resolution, i.e., number of points.
#############################################################

$r_num[0]=9;
$r_dr[0]=($r_max-$r_min)/$r_num[0];

print "number of resolution(r,t) : ", $num_dr, "\n";

#############################################################
#other resolutions, 
#each one is 2 time bigger than the previous one.
#############################################################

for ($i=1;$i<$num_dr;$i++)
{
    $r_num[$i]=$r_num[$i-1]*2;
    $r_dr[$i] =($r_max-$r_min)/$r_num[$i];
}


print "reslutions (r,dr): ";

for ($i=0;$i<$num_dr;$i++)
{
    print  $r_num[$i], ' ', $r_dr[$i], "; ";
    
}
print "\n";

#############################################################
# initialize mid location r.
#############################################################

for ($i=0;$i<$num_dr;$i++)
{
    for ($j=0; $j<$r_num[$i]; $j++)
    {
	$r[$j][$i] = $r_min+$j*$r_dr[$i]+ $r_dr[$i]/2.0;
	#print $r[$j][$i], "\n";
    }
}

#############################################################
# initialize number counter of each portion.
#############################################################

for ($i=0; $i<$num_dr;$i++)
{
    for ($j=0; $j<$r_num[$i];$j++)
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
	    if($data[7]==0)
	    {
		$x[$num_p] = $data[0];
		$y[$num_p] = $data[1];
		$num_p++;
	    }
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
# i.e., probability at each distance.
# number density of particls is num_p/V_eff
# periodic images have to be considered.
#############################################################

	$coeff = $num_density_r/2.0/$pi/$num_o;
	for($p=0;$p<$num_o; $p++)
	{  
	    for($q=0;$q<$num_t;$q++)
	    {
		$x_12 = $x_o[$p]-$x_t[$q];
		$y_12 = $y_o[$p]-$y_t[$q];
		$r_12 = sqrt($x_12**2+$y_12**2);
		
		#print "r_12: ", $r_12, "\n";
		
		for ($i=0;$i<$num_dr;$i++)
		{
		    
		    $r_idx=($r_12-$r_min)/$r_dr[$i];
		    
                    ###########################################
                    # exclude two particles too close or far.
                    ###########################################
		    if($r_idx<0 || $r_idx >= $r_num[$i])
		    {
			next;
		    }
		    
		    ############################################
		    ############################################
		    $num[$r_idx][$i] += ($coeff/$r[$r_idx][$i]/$r_dr[$i]);
		} # i < num_dr
	    } # q < num_t
	} # p< num_o
	
#############################################################	
# increase counter of files
#############################################################
	$num_file ++;
	
    }
    
}

print "number of files processed: ", $num_file, "\n";

if ( $num_file > 0 ) 
{
    for ($i=0; $i<$num_dr; $i++)
    {
	for ($j=0; $j<$r_num[$i]; $j++)
	{
	    $num[$j][$i] /= $num_file;
	    
	}
    }
}

for ($i=0; $i<$num_dr; $i++)
{
    
    $file_out_name = ">".$file_out_prefix."_".$gap."_".$r_dr[$i].".dat";
    open(OUT, $file_out_name);
    
    for ($j=0;$j<$r_num[$i];$j++)
    {
	print OUT $r[$j][$i],' ', $num[$j][$i], "\n";
    }
    
    
    print "writing output file : ", $file_out_name, "\n";
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

