#MPI 

N=2520
M=1920
process=( 4 9 16 25 36 49 64 81 100 )
loops=( 10 50 100 200 500 1000 2000 5000 8000 15000 60000 100000 )
size=( 1 2 4 )
for k in "${size[@]}";do
	for j in "${loops[@]}";do
		for i in "${process[@]}";do
			n=`echo $N*$k | bc`
			m=`echo $M*$k | bc`
			mpiexec -n $i -f machines build/conv -r $n -c $m -m $j 1> /dev/null 2>out
			echo Process: $i - Loops: $j - Size: $n x $m `cat out` >> MPI_Results
			echo ""
		done 
	done
done

rm out