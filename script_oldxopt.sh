source /opt/anaconda2/bin/activate root
source activate my_root

#echo "=====Varying MN quantity on original====="
#echo "-----10 MNs-----"
#for run in {1..5}
#do
#    python2 simulation.py 10 200
#done
#
#echo "-----50 MNs-----"
#for run in {1..5}
#do
#    python2 simulation.py 50 200
#done
#
#echo "-----100 MNs-----"
#for run in {1..5}
#do
#    python2 simulation.py 100 200
#done
#
#echo "-----400 MNs-----"
#for run in {1..5}
#do
#    python2 simulation.py 400 200
#done
#
#echo "=====Varying MN quantity on optimized====="
#echo "-----10 MNs-----"
#for run in {1..5}
#do
#    python2 simulation_opt.py 10 200
#done
#
#echo "-----50 MNs-----"
#for run in {1..5}
#do
#    python2 simulation_opt.py 50 200
#done
#
#echo "-----100 MNs-----"
#for run in {1..5}
#do
#    python2 simulation_opt.py 100 200
#done
#
#echo "-----400 MNs-----"
#for run in {1..5}
#do
#    python2 simulation_opt.py 400 200
#done

echo "=====Varying time on original====="
echo "-----200 ms-----"
for run in {1..5}
do
    python2 simulation.py 50 200
done

echo "-----500 ms-----"
for run in {1..5}
do
    python2 simulation.py 50 500
done

echo "-----1000 ms-----"
for run in {1..5}
do
    python2 simulation.py 50 1000
done

echo "=====Varying time on optimized====="
echo "-----200 ms-----"
for run in {1..5}
do
    python2 simulation_opt.py 50 200
done

echo "-----500 ms-----"
for run in {1..5}
do
    python2 simulation_opt.py 50 500
done

echo "-----1000 ms-----"
for run in {1..5}
do
    python2 simulation_opt.py 50 1000
done

source deactivate
