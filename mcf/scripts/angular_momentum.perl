open (fileinput, "mcf_particles.out");

$n = 0;
$t = 0.0;

$Lx=4.0;
$Ly=4.0;
$x0=$Lx/2.0;
$y0=$Ly/2.0;
$rho=1.0;
$m_t= $rho*$Lx*$Ly;

# total momentum inertia
$mmi_t = $m_t*($Lx+$Ly)/12.0;

#total momentum
$mom_t = 0.0;

while($line = <fileinput>)
{
    @data = split(' ', $line); 
    
    $x = $data[0];
    $y = $data[1];
    $vx = $data[2];
    $vy = $data[3];
    $m  = $data[5];
    
    #print $x, ' ', $y, ' ' , $vx, ' ', $vy, ' ', $m, "\n";
    
    $rx = ($x-$x0);
    $ry = ($y-$y0);
    $r2 = $rx**2 + $ry**2;
    $v2 = $vx**2 + $vy**2;
    
    $mom = sqrt($r2*$v2 - ($rx*$vx + $ry*$vy)**2);
    $mom *= $m;
    if ($rx*$vy-$ry*$vx > 0 )
    {
	$mom_t += $mom;
    }
    else
    {
	$mom_t -= $mom;
    }
}

print "angular momentum total     : ", $mom_t, "\n";
print "momentum of inertia total  : ", $mmi_t, "\n";
print "avarage angular velocity   : ", $mom_t/$mmi_t, "\n";

close(fileinput);

