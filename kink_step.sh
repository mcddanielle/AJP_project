#Modify DC current and run the MD simulation by modifying the input file Pcw0

DIR1=""
echo $DIR1

for current in $(seq 0.01 0.1 0.01)
do
    echo $current
    dir0=current_'printf "%1.2" $current'  #check number precision
    echo $dir0

    #-------------------
    if [ ! -d $dir0 ]
    then
        mkdir $dir0
    fi
    #------------------

    cd $dir0
    pwd

    restart=0

    if [ $restart -eq 0 ] && [ -f *.sh.o* ]
    then
	echo "Already run"
    else
	#make the new configuration file Pa0
	sed 's/current 0.2/current '"$current"'/' <../Pcw0> Pcw0
	cp ../restart_config .
	#run the physics MD code
	rm velocity_data/XV_data_t\=00*

	#run the simulation
	/home/danielle/Research/YDrive_AR/Source_Code/run_wstripe

	#extract velocity info into local file
	python measuring_kink_velocity.py > avg_kink_vel.txt

	#extract the first line of the file
	avg=`tail -n 1 avg_kink_vel.txt`

    fi

    cd ../

    #combine all info in a single file
    echo $current $avg >> avg_kink_vel.txt
done
	
exit
