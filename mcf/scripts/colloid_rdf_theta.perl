#############################################################
# Calculate radial distribution function of colloid pair
# according to relative angle.
#############################################################
use Math::Trig;

$dir_name=$ARGV[0];
$file_in_prefix="mcf_colloid";
$file_out_prefix=$ARGV[1];

opendir(DIR, $dir_name) || 
die ("can not open directory  :".$dir_name.", as it does not exist !\n");
print "processing directory   : ", $dir_name, "\n";

@file_names_temp = readdir(DIR);
@file_names_in=sort @file_names_temp;
closedir(DIR);
print "number of files inside : ", scalar(@file_names_in), "\n";


$pi=3.1415926;
$num_file  = 0;
$step_start=80000;
$step_start=$ARGV[2];
#$step_end=1000000;
$step_end  =99999999;
$step_end=$ARGV[3];

######################################
# Lx,Ly: box size.
# V    : volume 3D / area 2D.
# gap  : gap considered away from the walls.
# suppose x-direction periodic and
# y-direction wall.
######################################

$Lx=$ARGV[4];
$Ly=$ARGV[5];
$V=$Lx*$Ly;
$gap=$ARGV[6];
$R=1.0;
$V_eff=$Lx*($Ly-2.0*($gap-$R));


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
# minimum/maximum distrance
# minimum angle=0; maximum angle=pi.
#############################################################
$r_min=2.0;
$r_max=2.05;
$dr=$r_max-$r_min;
$t_min=0.0;
$t_max=$pi;

#############################################################
# number of different resolutions in angle.
#############################################################

$num_dt=5;

#############################################################
#the 0th resolution, i.e., number of points.
#############################################################


$t_num[0]=5;
$t_dt[0]=($t_max-$t_min)/$t_num[0];

print "number of resolution: ",  $num_dt,"\n";

#############################################################
#other resolutions
#############################################################

for ($i=1;$i<$num_dt;$i++)
{
    $t_num[$i]=$t_num[$i-1]*2;
    $t_dt[$i] =($t_max-$t_min)/$t_num[$i];
}

print "reslutions (t,dt): ";

for ($i=0;$i<$num_dt;$i++)
{
    print $t_num[$i], ' ',$t_dt[$i],"; "
}

print "\n";


#############################################################
# initialize mid location t.
#############################################################

for ($i=0;$i<$num_dt;$i++)
{
    for ($j=0; $j<$t_num[$i]; $j++)
    {
	$t[$j][$i] = $t_min+$j*$t_dt[$i]+ $t_dt[$i]/2.0;
	#print $t[$j][$i], "\n";
    }
}


#############################################################
# initialize number counter of each portion.
#############################################################

for ($i=0; $i<$num_dt;$i++)
{
    for ($j=0; $j<$t_num[$i];$j++)
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
print "starting file : ", $f_start, "\n";
print "ending file   : ", $f_end, "\n";

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
	
	#print "volume      : ", $V, "\n";
	#print "volume effective : ", $V_eff, "\n";
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
#############################################################
	$area  = $pi*($r_max**2-$r_min**2);
	$coeff = $num_o*$area*$num_density;
	
	#print "area        :", $area, "\n";
	#print "num_density :", $num_density, "\n";
	#print "coefficient :", $coeff, "\n";
	
	for($p=0; $p<$num_o; $p++)
	{
	    for($q=0; $q<$num_t; $q++)
	    {
		$x_12 = $x[$p]-$x[$q];
		$y_12 = $y[$p]-$y[$q];
		$r_12 = sqrt($x_12**2+$y_12**2);
		
		if ($r_12>$r_max || $r_12 < $r_min )
		{
		    next;
		}
		
		$t_12 = acos($x_12/$r_12);
		
		if ($y_12 < 0.0)
		{
		    $t_12 = $pi - $t_12;
		}
		
		#print "t_12: ", $t_12, "\n";
		
		for ($i=0; $i<$num_dt; $i++)
		{
		    $t_idx = ($t_12-$t_min)/$t_dt[$i];
		    $num[$t_idx][$i] += $t_num[$i]/$coeff;
		    
		} # i < num_dt
		
	    }# q < num_t
	    
	}# p < num_o
    
#############################################################	
# increase counter of files
#############################################################
	$num_file ++;
	
    } # if f is in consideration
    
}# folder

print "number of files processed: ", $num_file, "\n";

if ( $num_file > 0 ) 
{
    for ($i=0; $i<$num_dt; $i++)
    {
	for ($j=0; $j<$t_num[$i]; $j++)
	{
	    $num[$j][$i] /= $num_file;
	    
	}
    }
}

for ($i=0; $i<$num_dt; $i++)
{
    $file_out_name = ">".$file_out_prefix."_".$gap."_".$t_dt[$i].".dat";
    open(OUT, $file_out_name);
    
    for ($j=0;$j<$t_num[$i];$j++)
    {
	print OUT $t[$j][$i], ' ', $num[$j][$i], "\n";
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

