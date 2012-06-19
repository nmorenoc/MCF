###################################################
#Calculate average stress on the wall
###################################################

open (file_in, "mcf_boundary.dat");


###################################################
# indicate if this script is used as batch processing
# if so, only data are output, other information
# are omitted.
###################################################
$batch_process=1;


$x_min=0;
$x_max=5.e5;

#####################
# time period of interest.
#####################
$t_start = 40.0;
$t_end   = 1000.0;

$t_start = $ARGV[0];
$t_end   = $ARGV[1];

#print "t_start:", $t_start, "\n";
#print "t_end  :", $t_end, "\n";


$t = 0.0;
$n = 0;

$down_x = 0.0;
$down_y = 0.0;
$up_x = 0.0;
$up_y = 0.0;

$down_x2 = 0.0;
$down_y2 = 0.0;
$up_x2 = 0.0;
$up_y2 = 0.0;


while($line = <file_in>)
{
    @data = split(' ', $line); 

    $t = $data[1];
    if( $t<=$t_end )
    {
	if( $t>=$t_start )
	{
	    if (abs($data[2])< $x_min || abs($data[2])> $x_max)
	    {
		next;
	    }
	    if (abs($data[4])< $x_min || abs($data[4])> $x_max)
	    {
		next;
	    }
	    $n = $n+1.0;
	    $down_x = $down_x + $data[2];
	    $down_x2 = $down_x2 + $data[2]**2;
	    $down_y = $down_y + $data[3];
	    $down_y2 = $down_y2 + $data[3]**2;
	    $up_x = $up_x + $data[4];
	    $up_x2 = $up_x2 + $data[4]**2;
	    $up_y = $up_y + $data[5];
	    $up_y2 = $up_y2 + $data[5]**2;	    
	}   
    }
    else
    {
	last;
    }
}

########################################
#n=0: indicate no data for processing.
########################################
if ( $n==0 )
{
    print "", 0.0, ' ', 0, ' ', 0, ' ';
    print "", 0.0, ' ', 0, ' ', 0, ' ';
    print "", $n, "\n";
}
else
{
    $down_x = $down_x/$n;
    $down_x2 = $down_x2/$n;
    $down_y = $down_y/$n;
    $down_y2 = $down_y2/$n;
    $up_x = $up_x/$n;
    $up_x2 = $up_x2/$n;
    $up_y = $up_y/$n;
    $up_y2 = $up_y2/$n;
    
    $down_sx = sqrt($down_x2-$down_x**2);
    $down_sy = sqrt($down_y2-$down_y**2);
    $up_sx=sqrt($up_x2-$up_x**2);
    $up_sy=sqrt($up_y2-$up_y**2);
    $down_ex = $down_sx/sqrt($n);
    $down_ey = $down_sy/sqrt($n);
    $up_ex = $up_sx /sqrt($n);
    $up_ey = $up_sy/sqrt($n);
    
    if ($batch_process) 
    {
	print "", ($down_x - $up_x)/2.0, ' ' , ($down_sx+$up_sx)/2.0, ' ' , ($down_ex+$up_ex)/2.0, ' ';
	print "", ($up_y - $down_y)/2.0, ' ' , ($up_sy+$down_sy)/2.0, ' ' , ($up_ey+$down_ey)/2.0, ' ';
	print "", $n, "\n";
    }
    else
    {
	print "n      : ", $n, "\n";
	print "stress, standard deviation, error\n";
	print "down_x,vari,error : ", $down_x,' ',$down_sx, ' ' , $down_ex, "\n";
	print "down_y,vari,error : ", $down_y,' ',$down_sy , ' ', $down_ey,"\n";
	print "up_x,vari,error   : ", $up_x, ' ',$up_sx,' ' , $up_ex,"\n";
	print "up_y,vari,error   : ", $up_y, ' ',$up_sy,' ' , $up_ey,"\n";
	print "Average :\n";
	print "x, vari,error     : ";
	print "", ($down_x - $up_x)/2.0, ' ' , ($down_sx+$up_sx)/2.0, ' ' , ($down_ex+$up_ex)/2.0, "\n";
	print "y, vari,error     : ";
	print "", ($up_y - $down_y)/2.0, ' ' , ($up_sy+$down_sy)/2.0, ' ' , ($up_ey+$down_ey)/2.0, "\n";
	print "", $n, "\n";
    }
} # n=0
close(file_in);
    
    
