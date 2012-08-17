#############################################################
# Build up punto animation file with different colors
# for colloid and orientation marker of colloid,
# and wall particles.
# NOTE:
# We assume that x direction is possibly periodic and
# fluid particles are removed.
#############################################################


use Math::Trig;

$R =$ARGV[0];
$Lx=$ARGV[1];
$Ly=$ARGV[2];
$Nc=$ARGV[3];
$folder_in=$ARGV[4];
print "R: ", $R, "\n";
print "Lx: ", $Lx, "\n";
print "Ly: ", $Ly, "\n";
print "Nc: ", $Nc, "\n";

#############################################################
#Open a directory and read in file names.
#############################################################

opendir (dir_in, $folder_in) || die("can not open directory");
print "Processing directory: ", $folder_in, "\n";

@names_in=readdir(dir_in);
#print scalar @names_in, "\n";

@names_in_sorted=sort @names_in;

#############################################################
#open output file.
#############################################################

$name_out=">punto_particles_marker_no_fluid.dat";
open file_out, $name_out;


#############################################################
#pid_index: the column of particle ID.
#sid_index: the column of species ID.
#
#color_* to indicate colors of
#solvent, wall, colloid, and marker of colloid.
#############################################################

$x_index=0;
$y_index=1;
$pid_index=6;
$sid_index=7;


$color_s=0;
$color_w=-1;
$color_c_min=1;
$color_c_max=1.8;
$color_m=2;

#############################################################
#We use first input file to make a hash table
#(key=pid, value=sid/color)
#which uses particle id(pid) to 
#determine color of each particle.
#
#Rotating x+ axis around (0,0) between a_min and a_max
#to mark a certain region of a colloid with a special color.
#This is useful for a rotating colloid, as it tells
#the orientation of each colloid during an animation.
#############################################################

$hash_made=0;
$a_min=-0.5;
$a_max=0.5;

#############################################################
#Since we have Nc colloids, reset its center position and
#number of boundary particles to be zero.
#We will cacluate center of each one
#by averaging all the boundary particels's positions later.
#############################################################

for ($i=1;$i<=$Nc;$i++)
{
    $x_col[$i]=0;
    $y_col[$i]=0;
    $Nb[$i]=0;
}

#############################################################
#Count how many boundary particles at minimum and maximum
#end of the box for each colloid, to determin the colloid
#cross periodic boundary.
#############################################################
for ($i=1;$i<=$Nc;$i++)
{
    $xc_min[$i]=0;
    $xc_max[$i]=0;
    $cross[$i]=0;
}


#############################################################
#read the first file in this folder.
#use the first file to make hash table.
#############################################################
	
foreach $f (@names_in_sorted)
{
    #########################################################
    # Select the first particles file to read.
    #########################################################
    
    if ( ($f =~ "particles") && ( $f =~ ".out" ) )
    {
	open file_in, $f;
	
	#####################################################
	#read in the frist file first time to determine
	#if the colloid cross a periodic boundary
	#####################################################
	
	while($line=<file_in>)
	{
	    @data=split(' ', $line);
	    
	    $sid=$data[$sid_index];
	    
	    if ( $sid>0 )
	    {
		if ( $data[$x_index] < 2.0*$R ) 
		{
		    $xc_min[$sid]+=1;
		}elsif ( $data[$x_index] > $Lx-2.0*$R ) 
		{
		    $xc_max[$sid]+=1;
		}
		
		$Nb[$sid]++;
		
		if($sid>$Nc)
		{
		    $Nc=$sid;
		}
	    }
	}#while (line=<file_in>)
	
	for ($i=1;$i<=$Nc;$i++)
	{
	    if ( $xc_min[$i] > 0 && $xc_max[$i] > 0 ) 
	    {
		$cross[$i]=1;
		#print "i, cross[i]", $i, ' ' , $cross[$i], "\n";
	    }
	}
	
	close(file_in);
	    
	#####################################################
	#read in the frist file second time 
	#to calculate colloid center
	#####################################################
	
	open file_in, $f;
	
	while($line=<file_in>)
	{
	    @data=split(' ', $line);
	    
	    $sid=$data[$sid_index];
	    
	    #################################################
	    #for each colloidal boundary particle
	    #if the colloid crossing periodic boundary
	    #we shift its x center to the middle of the box,
	    #which is necessary to calculate its oriention
	    #by marker colloidal boundary particles.
	    #################################################
	    if ( $sid>0 )
	    {
		if ( $cross[$sid] == 1 && $data[$x_index] < 2.0*$R ) 
		{
		    $data[$x_index]+=$Lx/2.0;
		    
		}elsif ( $cross[$sid] == 1 && $data[$x_index] > $Lx-2.0*$R ) 
		{
		    $data[$x_index]-=$Lx/2.0;
		}
		
		$x_col[$sid]+=$data[$x_index];
		$y_col[$sid]+=$data[$y_index];
		
	    }
	}#while (line=<file_in>)
	
	#####################################################
	#estimate center of each colloid, note that the one
	#crossing periodic boundary, its center is shifted.
	#####################################################
	
	for ($i=1;$i<=$Nc;$i++)
	{
	    $x_col[$i]/=$Nb[$i];
	    $y_col[$i]/=$Nb[$i];
	    #print $x_col[$i], ' ',$y_col[$i], "\n";
	}
	
	close(file_in);
	
	#####################################################
	#read the third times this file to make hash table
	#mapping particle id to color.
	#####################################################
	
	open file_in, $f;
	
	while($line=<file_in>)
	{
	    @data=split(' ', $line);
		
	    $pid=$data[$pid_index];
	    $sid=$data[$sid_index];
	    
	    #################################################
	    #for a colloidal boundary particle, 
	    #scale its color according to min, max color
	    #and number of total colloids.
	    #Also determine if it is for indicating
	    #orientation.
	    #################################################
	    if ( $sid>0 )
	    {
		if ($Nc>1)
		{
		    $color_c = $color_c_min+
			($sid-1)*($color_c_max-$color_c_min)/($Nc-1);
		    
		}else
		{
		    $color_c = $color_c_min;
		}
		#############################################
		#shift the boundary particle, if the 
		#corresponding colloid cross periodic boundary.
		#############################################
		
		if ( $cross[$sid] ==1 && $data[$x_index] < 2.0*$R ) 
		{
		    $data[$x_index]+=$Lx/2.0;

		}elsif ( $cross[$sid] ==1 && $data[$x_index] > $Lx-2.0*$R ) 
		{
		    $data[$x_index]-=$Lx/2.0;
		}
		
		#############################################
		#consider the right half the colloid.
		#############################################

		if ( $data[$x_index]>$x_col[$sid] )
		{
		    $angle=
			atan(($data[$y_index]-$y_col[$sid])/
			     ($data[$x_index]-$x_col[$sid]));
		}
		else 
		{
		    $angle = -10.0;
		}
		
		if($angle >= $a_min && $angle <= $a_max )
		{
		    $hash_ps{$pid}=$color_m;
		}
		else
		{
		    $hash_ps{$pid}=$color_c;
		}
		
	    }
	    #################################################
	    #for wall boundary particle, simple use
	    #assgined color.
	    #################################################
	    elsif( $sid < 0 )
	    {
		$hash_ps{$pid}=$color_w;
	    }
	    
	}#while
	
	close(file_in);
	
	#################################################
	#hash table has been made.
	#################################################
	last;
    }    
}

#############################################################
#once hash tale for color is built,
#read all files in this folder.
#############################################################

foreach $f (@names_in_sorted)
{
    if ( ($f =~ "particles") && ( $f =~ ".out" ) )
    {
	print $f, ">", $name_out, "\n";
	
	open file_in, $f;
	while($line=<file_in>)
	{
	    @data=split(' ', $line);
	    $pid=$data[$pid_index];
	    $sid=$data[$sid_index];
	    
	    #################################################
	    #print x, y, vx, vy, and color.
	    #################################################
	    
	    if ( $sid != 0 ) 
	    {
		for($i=0;$i<4;$i++)
		{
		    print file_out $data[$i], ' ';
		}
		print file_out $hash_ps{$pid}, "\n";
	    }
	}
	print file_out "\n";
	
	
	close(file_in);
    }
}


closedir(dir_in);
close(file_out);
