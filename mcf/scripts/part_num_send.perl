open (fileinput, "part_num_send.dat");

$n = 0;

@send[32]=0;
@recv[32]=0;
@sendn[32]=0;
@recvn[32]=0;

while($line = <fileinput>)
{
    @data = split(' ', $line); 
    
  if ($data[1]>=0) 
    {

   	$send[$data[0]]  += $data[2];
	$sendn[$data[0]] += 1;
	
	$recv[$data[1]]  += $data[2];
       	$recvn[$data[1]] += 1;
	
	if($data[1]==2) 
	{
	    print $data[1], ' receive from ', $data[0], ' ' , $data[2], "\n";
	}
	if($data[0]==10) 
	{
	    print $data[0], ' send to ', $data[1], ' ', $data[2], "\n";
	}
    }  
    
}

for($i=0;$i<32;$i++)
{
#    print 'rank, send, recv  : ', $i, ' ', $send[$i], ' ', $recv[$i], "\n";
    #print 'rank, sendn,recvn ; ', $i, ' ', $sendn[$i], ' ', $recvn[$i], "\n";
    print  $i, ' ', $sendn[$i], ' ', $recvn[$i], "\n";

}

close(fileinput);

