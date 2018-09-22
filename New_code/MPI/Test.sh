#MPI 
M=2520
N=1920
make mpi
process=( 4 9 16 25 36 49 64 81 100 )
loops=( 100 500 1000 2000 5000 10000 15000 60000 100000 1000000 )
size=( 0.25 0.5 1 2 4 )
for k in "${size[@]}";do
	for j in "${loops[@]}";do
		for i in "${process[@]}";do
			row=`echo $N*$k/1 | bc`
			col=`echo $M*$k/1 | bc`
			mpiexec -n $i -f machines build/conv -r $row -c $col -m $j -f  $1 2> /dev/null 1>out
			echo Process: $i - Loops: $j - Size: $row x $col `cat out` >> MPI_Results
			echo "" >> MPI_Results
		done
		echo "" >> MPI_Results
		echo "" >> MPI_Results
		echo Finished loops: $j 
	done
done

for k in "${size[@]}";do
	for j in "${loops[@]}";do
		for i in "${process[@]}";do
			row=`echo $N*$k/1 | bc`
			col=`echo $M*$k/1 | bc`
			mpiexec -n $i -f machines build/conv -r $row -c $col -m $j -f  $2 -rgb 2> /dev/null 1>out
			echo Process: $i - Loops: $j - Size: $row x $col `cat out` >> MPI_RGB_Results
			echo "" >> MPI_RGB_Results
		done
		echo "" >> MPI_RGB_Results
		echo "" >> MPI_RGB_Results
		echo Finished loops: $j 
	done
done

make openmp
process=( 4 9 16 25 )
loops=( 100 500 1000 2000 5000 10000 15000 60000 100000 1000000 )
size=( 0.25 0.5 1 2 4 )
threads=( 2 3 4 5 6 7 )
for k in "${size[@]}";do
	for j in "${loops[@]}";do
		for i in "${process[@]}";do
			for t in "${threads[@]}";do
				row=`echo $N*$k/1 | bc`
				col=`echo $M*$k/1 | bc`
				mpiexec -n $i -f machines build/conv -r $row -c $col -m $j -f  $1 -t $t 2> /dev/null 1>out
				echo Process: $i - Threads: $t - Loops: $j - Size: $row x $col `cat out` >> Hybrid_Results
				echo "" >> Hybrid_Results
			done
		done
		echo "" >> Hybrid_Results
		echo "" >> Hybrid_Results
		echo Finished loops: $j 
	done
done

for k in "${size[@]}";do
	for j in "${loops[@]}";do
		for i in "${process[@]}";do
			for t in "${threads[@]}";do
				row=`echo $N*$k/1 | bc`
				col=`echo $M*$k/1 | bc`
				mpiexec -n $i -f machines build/conv -r $row -c $col -m $j -f  $2 -t $t -rgb 2> /dev/null 1>out
				echo Process: $i - Threads: $t - Loops: $j - Size: $row x $col `cat out` >> Hybrid_RGB_Results
				echo "" >> Hybrid_RGB_Results
			done
		done
		echo "" >> Hybrid_RGB_Results
		echo "" >> Hybrid_RGB_Results
		echo Finished loops: $j 
	done
done

rm out