
###################################################
# Mark species ID plus one extra line of boundary
# particles to indicate rotation of 
# a colloidal particle.
# NOTE:
# all colloidal particle should be inside compuational
# domain, since we want to use all boundary particles
# to estimate their centers.
###################################################

#####################
#Open a directory.
#read in file names.
#####################

use Math::Trig;
$folder_in=$ARGV[0];
opendir (dir_in, $folder_in) || die("can not open directory");
#print "Processing directory: ", $folder_in, "\n";

@names_in=readdir(dir_in);
#print scalar @names_in, "\n";
#print  @names_in, "\n";

@names_in_sorted=sort @names_in;

#####################
#open a file for output.
#####################

$name_out=">punto_particles_marker.dat";
open file_out, $name_out;



#####################
#sid/color to make.
# solvent, colloid
# and maker.
#####################
$sid_w=-1;
$sid_s=0;
$sid_c_min=1;
$sid_c_max=1.2;
$sid_m=1;

$pid_index=8;
$sid_index=9;


#####################
#we use first input file
#to make a hash table
#(key=pid, value=sid/color)
#which uses particle id
#to tell sid/color of
#each boundary particle.
#
#boundary line angle between
# a_min and a_max
# are marked with
# special color;
#####################
$hash_made=0;
#%hash_ps;
$a_min=-0.5;
$a_max=0.5;

#####################
#Suppose we have Np colloids
#cacluate center of each one
#by averaging all the boundary
#particels's positions.
#####################
$Nc=200;
for ($i=1;$i<=$Nc;$i++)
{
    $x_col[$i]=0;
    $y_col[$i]=0;
    $z_col[$i]=0;
    $Nb[$i]=0;
}
$Nc=0;

#####################
#read each file in this folder.
#####################

foreach $f (@names_in_sorted)
{
    if ( ($f =~"particles") && ( $f =~ ".out" ) )
    {
	print $f, ">", $name_out, "\n";
	#open file
	open file_in, $f;
	
	#use the first file to make hash table.
	if ($hash_made==0)
	{
	    #read the file once to calcualte the center.
	    while($line=<file_in>)
	    {
		@data=split(' ', $line);
		
		$sid=$data[$sid_index];
		#collodal boundary particle
		if ($sid>0)
		{
		    if($sid>$Nc)
		    {
			$Nc=$sid;
		    }
		    $x_col[$sid]+=$data[0];
		    $y_col[$sid]+=$data[1];
		    $z_col[$sid]+=$data[2];
		    $Nb[$sid]++;		    
		}
	    }#while
	    
	    #estimate center of each colloid.
	    for ($i=1;$i<=$Nc;$i++)
	    {
		$x_col[$i]/=$Nb[$i];
		$y_col[$i]/=$Nb[$i];
		$z_col[$i]/=$Nb[$i];
		#print $x_col[$i], ' ',$y_col[$i], ' ',$z_col[$i], "\n";
	    }
	    
	    close(file_in);
	    open file_in, $f;
	    #read the again to make hash table
	    while($line=<file_in>)
	    {
		@data=split(' ', $line);
		$pid=$data[$pid_index];
		$sid=$data[$sid_index];
		if ($Nc>1)
		{
		    $sid_c = $sid_c_min+
			($sid-1)*($sid_c_max-$sid_c_min)/($Nc-1);
		    #if ($sid==1)
		    #{$sid_c=0;}
		}else
		{
		    $sid_c = $sid_c_min;
		}
		if ($sid==0)
		{
		    $hash_ps{$pid}=$sid_s;
		}elsif($sid>0)
		{
		    $angle=$a_min-1;
		    
		    if ($data[0]>$x_col[$sid])
		    {
			$angle=
			    atan(($data[1]-$y_col[$sid])/
				 ($data[0]-$x_col[$sid]));
		    }
		    
		    if($angle >= $a_min && $angle <= $a_max )
		    {
			$hash_ps{$pid}=$sid_c+1;
		    }
		    else
		    {
			$hash_ps{$pid}=$sid_c;
		    }
		    
			
		}else
		{
		    $hash_ps{$pid}=$sid_w;
		}
		for($i=0;$i<6;$i++)
		{
		    print file_out $data[$i], ' ';
		}
		print file_out $hash_ps{$pid}, "\n";
	    }#while
	    print file_out "\n";	    
	    $hash_made=1;
	}#if (hash_made==0)
	else
	{
	    while($line=<file_in>)
	    {
		@data=split(' ', $line);
		$pid=$data[$pid_index];
		$sid=$data[$sid_index];
		for($i=0;$i<6;$i++)
		{
		    print file_out $data[$i], ' ';
		}
		print file_out $hash_ps{$pid}, "\n";
	    }
	    print file_out "\n";
		
	}

	close (file_in);

    }
}

closedir(dir_in);
close(file_out);
