for i in {44..50}
do
    ipython zp_rotate.py $i
    python simulate_zp_after_tilt.py $i
done
mkdir results
mv *.npy results/
