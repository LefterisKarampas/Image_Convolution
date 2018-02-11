rm machines
for i in `seq 1 30`;
        do
        	if [ $i -lt 10 ]
        		then
        		ssh linux0$i date
        		if [ $? -eq 0 ]
					then
						echo "linux0$i:$1" >> machines
					fi
        	else
        		ssh linux$i date
        		if [ $? -eq 0 ]
				then
					echo "linux$i:$1" >> machines
				fi
        	fi
        done    
