$n = 0;

for ($i=6; $i<=7; $i ++ )
{
    for ($j=6; $j<=7 ; $j++)
    {
	for ($k= 12; $k<=13; $k++)
	{
	    $n ++;
	    $num_part = $i * $j * $k;
	    print 'n        :', $n, "\n";
	    print 'num_part : ', $num_part, "\n";
		
	}
    }
	
}
