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


    #extract velocity info into local file
    python3 ~/pymodules/animation_code/measuring_kink_velocity.py > avg_kink_vel.txt

    #extract the first line of the file
    avg=`tail -n 1 avg_kink_vel.txt`
	
    cd ../
    pwd

    #combine all info in a single file
    echo $current $avg >> avg_kink_vel.txt
done
	

