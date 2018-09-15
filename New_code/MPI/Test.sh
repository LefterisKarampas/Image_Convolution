#MPI 
M=2520
N=1920
process=( 4 9 16 25 36 49 64 81 100 )
## loops=( 100 500 1000 2000 5000 ) ## 10000 15000 60000 100000)
loops=( 100000 )
size=( 1 )
for k in "${size[@]}";do
	for j in "${loops[@]}";do
		for i in "${process[@]}";do
			row=`echo $N*$k | bc`
			col=`echo $M*$k | bc`
			mpiexec -n $i -f machines build/conv -r $row -c $col -m $j -f  $1 -o $2 2> /dev/null 1>out
			echo Process: $i - Loops: $j - Size: $row x $col `cat out` >> MPI_Results
			echo "" >> MPI_Results
		done
		echo "" >> MPI_Results
		echo "" >> MPI_Results
		echo Finished loops: $j 
	done
done

rm out