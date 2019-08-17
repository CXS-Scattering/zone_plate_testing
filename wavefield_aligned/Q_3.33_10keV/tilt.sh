for i in {33..50}
do
    ipython zp_rotate.py $i
    python simulate_zp_with_tilt.py $i
done
