#Modify DC current and run the MD simulation by modifying the input file Pcw0

DIR1=""
echo $DIR1

for current in $(seq 0.001 0.001 0.1)
do
    echo $current
    string_current=`printf "%1.3f" $current`
    dir0=current_$string_current  #check number precision
    echo $dir0

    #-------------------
    if [ ! -d $dir0 ]
    then
        mkdir $dir0
    fi
    #------------------

    cd $dir0
    pwd


    #make the new configuration file Pa0
    sed 's/current 0.2/current '"$current"'/' <../Pcw0> Pcw0

    #the high DC currents (current > 0.5) don't need to run
    #for so long we could include another sed call to reduce the total simulation time
    #
    
    
    cp ../restart_config .
    #run the physics MD code
    rm velocity_data/XV_data_t\=00*

    #run the simulation
    #/home/danielle/Research/YDrive_AR/Source_Code/substrate_run_wstripe

    #extract velocity info into local file
    python3 ~/pymodules/animation_code/measuring_kink_velocity.py > avg_kink_vel.txt

    #extract the first line of the file
    #avg=`tail -n 1 avg_kink_vel.txt`

    #python ~/pymodules/animation_code/channel_colloid_movie_maker.py 
	
    cd ../
    pwd

    #combine all info in a single file
    echo $current $avg >> avg_kink_vel.txt
done
	

