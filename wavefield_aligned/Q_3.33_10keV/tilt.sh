for i in {1..50}
do
    ipython zp_rotate.py $i
    python simulate_zp_with_tilt.py $i
done
mkdir results
mv *.npy results/
