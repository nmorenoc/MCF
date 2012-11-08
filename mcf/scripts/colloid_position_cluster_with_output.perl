#############################################################
# Calculate size distribution of colloid-clusters by positions
# averaged by time.
# It is the same as colloid_position_cluster.perl
# However, there may be two more output from this perl script.
#          1:) output cluster ID back to colloid files.
#          2:) output size distribution of colloid-clusters
#              at any time instant.
#############################################################
use Math::Trig;


#############################################################
#freq_cluster_output_ID: how frequent is the .out with cluster
#                        ID is output. 
#                        0: no output at all.
#############################################################
$freq_cluster_ID_output = 1;
$dir_name_cluster_ID_output="cluster_ID_output/";

#############################################################
#freq_probability_cluster_output: how frequent to output
#                                 probability of cluster-size
#                                 at certain time step.
#############################################################
$freq_cluster_probability_output = 1;
$dir_name_cluster_probability_output="cluster_probability_output/";

$dir_name=$ARGV[0];
$file_in_prefix="mcf_colloid";
$file_out=$ARGV[1];

opendir(DIR, $dir_name) || 
die   ("can not open directory   :".$dir_name.", as it does not exist !\n");
print "processing directory      : ", $dir_name, "\n";

@file_names_temp = readdir(DIR);
@file_names_in=sort @file_names_temp;
closedir(DIR);
print "number of files inside    : ", scalar(@file_names_in), "\n";

$pi=3.1415926;
$num_file  = 0;
$step_start=80000;
$step_start=$ARGV[2];
$step_end  =99999999;
$step_end=$ARGV[3];

#############################################################
# Lx,Ly: box size.
# V    : area 2D.
# gap  : gap considered away from the walls.
# suppose x-direction periodic and y-direction wall.
#############################################################

$Lx=$ARGV[4];
$Ly=$ARGV[5];
$V=$Lx*$Ly;
$N=192;
$N=$ARGV[6];
$gap=$ARGV[7];
$R=1.0;
$V_eff=$Lx*($Ly-2.0*($gap-$R));
print "Lx, Ly, N, gap, R: ", $Lx, '; ', $Ly, '; ', $N, '; ', $gap, '; ', $R, ".\n";

#############################################################
# Surface distance threshold for defining cluster.
#############################################################
$cluster_s = 0.05;
$cluster_s = $ARGV[8];
print "cluster_s: ", $cluster_s, "\n";

#############################################################
# two dimensional matrix to indicate near-by relation.
# 1: related
# 0: unrelated.
#############################################################

for ( $i=1; $i<=$N; $i++ )
{
    for ( $j=1; $j<=$N; $j++ )
    {
	$near_by[$i][$j] = 0;
    }
}

#############################################################
# for recording the averaged number of different sized clusters.
#
# $num_cluster_size[$i] stores how many clusters with size $i.
# 
#############################################################

for ( $i=1; $i<=$N; $i++ )
{
    $num_cluster_size[$i]=0.0;
}

#############################################################
#num_p: number of colloids in the file.
#num_o: number of colloids considered,  others are exclused,
#       such as the ones near walls.
#num_o_average: averaged number of colloids considered 
#               over time.
#############################################################

$num_p=0;
$num_o=0;
$num_o_average=0;

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
# reset counter $num_p to zero. 
# read each colloid position.
#############################################################
	
	open (IN, $file_name);	
	
	$num_p=0;
	
	while ($line = <IN>)
	{
	    @data = split(' ', $line);
	    
	    $num_p++;
	    $x[$num_p] = $data[0];
	    $y[$num_p] = $data[1];
	
	}
	close(IN);
	#print "num_p : ", $num_p, "\n";
	
#############################################################
# get number denisty.
#############################################################

	$num_density=$num_p/$V;
	$num_density_r=$V/$num_p;
	#print "num density : ", $num_density, "\n";
	
#############################################################
# reset counter $num_o to zero. 
# record the ones, which are used for computation.
#############################################################

	$num_o=0;

	for ($p=1; $p<=$num_p; $p++)
	{
	    if ( ( $y[$p]>$gap ) && ( $y[$p]< $Ly-$gap) )
	    {
		$num_o++;
		$x_o[$num_o] = $x[$p];
		$y_o[$num_o] = $y[$p];
		
	    }
	   
	}
	#print "num_o : ", $num_o, "\n";
	$num_o_average += $num_o;
	
#############################################################
# Decide whether a pair is considered as near-by by
# checking the surface distance.
#############################################################
	
	for ( $i=1; $i<$num_o; $i++ )
	{
	    for ( $j=$i+1; $j<=$num_o; $j++ )
	    {
		
#############################################################
# reset near-by to zero, no relation for each pair initially.
#############################################################
		$near_by[$i][$j] = 0;
		
#############################################################
# check the surface distance of the pair to 
# decide if they are near-by.
#############################################################
		$r_ij = sqrt(($x_o[$i]-$x_o[$j])**2 + ($y_o[$i]-$y_o[$j])**2);

		if ( $r_ij <= 2.0*$R + $cluster_s )
		{
		    $near_by[$i][$j] = 1;
		}
		else
		{
#############################################################
# check the surface distance of the pair by considering 
# periodic image to decide if they are near-by.
#############################################################

		    $r_ij = sqrt(($x_o[$i]-$x_o[$j]-$Lx)**2 + ($y_o[$i]-$y_o[$j])**2);
		    
		    if ( $r_ij <= 2.0*$R + $cluster_s )
		    {
			$near_by[$i][$j] = 1;
		    }
		    else
		    {
#############################################################
# check the surface distance of the pair by considering 
# periodic image to decide if they are near-by.
#############################################################
			
			$r_ij = sqrt(($x_o[$i]-$x_o[$j]+$Lx)**2 + ($y_o[$i]-$y_o[$j])**2);
			
			if ( $r_ij <= 2.0*$R + $cluster_s )
			{
			    $near_by[$i][$j] = 1;					    
			}
		    }
		}
	    } # j
	}# i

#############################################################
# Check two dimensional near_by matrix to find cluster.
# 1: reset each particle as un-assigned to any cluster.
#############################################################
	
	for ( $i=1; $i<=$num_o; $i++ )
	{
	    $assigned[$i] = 0;
	}
	
#############################################################
#2: assign everty particle in each cluster with 
#   the minimum particle ID within this cluster, i.e.,
#   $mark is given the minimum ID to mark each particle
#   within this cluster.
#############################################################

	for ( $i=1; $i<=$num_o; $i++ )
	{
	    if ( $assigned[$i] > 0 )
	    {
		$mark = $assigned[$i];
	    }
	    else
	    {
		$assigned[$i]=$i;
		$mark = $i;
	    }
		
	    for ( $j=$i+1; $j<=$num_o; $j++ )
	    {
		if ( $near_by[$i][$j] >0 ) 
		{
		    $assigned[$j] = $mark;
		}
	    }
	}
	
#############################################################	
# 2.5 write certain snapshot with cluster ID for animation.
#############################################################	
	if ($freq_cluster_ID_output !=0 && $num_file % $freq_cluster_ID_output == 0)
	{
	    $file_name_cluster_ID_output = ">". $dir_name. $dir_name_cluster_ID_output . $f;
	    #print $file_name_cluster_output_ID, "\n";
	    open (OUTPUT_ID, $file_name_cluster_ID_output);
	    for ( $j=1; $j<=$num_o; $j++ )
	    {
		print OUTPUT_ID $x_o[$j], ' ', $y_o[$j], ' ', $assigned[$j], "\n";
	    }
	    close(OUTPUT_ID);	     
	}
#############################################################	
#3 Counter number of particles in each cluster.
#############################################################
	
	for ( $i=1; $i<=$num_o; $i++)
	{
	    $num_particle_cluster[$i] = 0;
	}
	for ( $i=1; $i<=$num_o; $i++)
	{
	    $num_particle_cluster[$assigned[$i]]++;
	}

	#print "cluster ID, number of particles: \n";
	#for ($i=1;$i<=$num_o;$i++)
	#{
	#    print $i, ': ', $num_particle_cluster[$i], "\n";
	#}
	#last;
	
#############################################################
# for recording the number of different sized clusters at
# current time step.
#
# $num_cluster_size_current[$i] stores how many clusters with size $i.
# 
#############################################################

	for ( $i=1; $i<=$num_o; $i++ )
	{
	    $num_cluster_size_current[$i]=0.0;
	}
	
	for ( $i=1; $i<=$num_o; $i++)
	{
	    $num_cluster_size_current[$num_particle_cluster[$i]] ++;
	    $num_cluster_size[$num_particle_cluster[$i]] ++;
	}

#############################################################
# 3.5 output probability of cluster size at certain time steps.
#############################################################
	if ($freq_cluster_probability_output !=0 && $num_file % $freq_cluster_probability_output == 0)
	{
	    $file_out_pre=substr($file_out,0,4);
	    $file_number = substr($f, 11,8);
	    $file_name_cluster_probability_output = ">". $dir_name. $dir_name_cluster_probability_output . $file_out_pre . '_' .$file_number. ".dat";
	    #print $file_name_probability_cluster_output, "\n";
	    open (OUTPUT_PROB, $file_name_cluster_probability_output);
	    for ( $i=1; $i<=$num_o; $i++ )
	    {
		print OUTPUT_PROB $i, ' ', $num_cluster_size_current[$i]/$num_o, "\n";
	    }
	    close(OUTPUT_ID);	     
	    
	}

#############################################################	
# increase counter of files
#############################################################
	$num_file ++;

	#if ($num_file > 1 )
	#{ last;}
	
    } # if file is in consideration.
    
}# folder

print "total number of files processed: ", $num_file, "\n";
if ( $num_file < 1 )
{
    print "number of files is zero !\n";
    exit;
}

#############################################################
# 
#############################################################
$num_o_average /= $num_file;
print "averaged number of particles used in each file: ", $num_o_average, "\n";

#############################################################
# Output file name.
#############################################################
$file_out=">". $file_out. ".dat";

#open(OUT, $file_out);
#print OUT "#cluster size, number of this sized cluster \n";

$check=0;

for ( $i=1; $i<=$N; $i++ )
{
#############################################################
# First averaged over number of files processed.
#############################################################
    $num_cluster_size[$i]/=$num_file;
    
#############################################################
# Check the if averaged total number of particles conserved.
#############################################################

    $check+= $i*$num_cluster_size[$i];    

#############################################################
# Check the probability of particles in each-sized cluster. 
#############################################################
    
    $num_cluster_size[$i] =  $num_cluster_size[$i] * $i / $num_o_average;
    
    #print OUT $i, ' ',$num_cluster_size[$i], "\n";
}

print "check sum of probability: ", $check, "\n";
#close(OUT);

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

