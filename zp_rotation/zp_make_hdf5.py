from mpi4py import MPI
import h5py
import numpy as np
import os

'''
make_zp_from_rings : make a zone plate from the rings which were created earlier.

Inputs : n - number of rings, grid_size

Outputs : a numpy array containing the zone plate
'''
def make_zp_from_rings(n, grid_size):
    pwd = os.getcwd()
    os.chdir('/home/sajid/packages/zone_plate_testing/rings')
    zp = np.zeros((grid_size, grid_size))
    for i in range(n):
        if i % 2 == 1:
            ring_ = np.load('ring_' + str(i) + '.npy')
            ring_ = tuple((ring_[0], ring_[1]))
            zp[ring_] = 1
    os.chdir(pwd)
    return zp


if __name__ == "__main__":
    # The process ID (integer 0-3 for 4-process run)
	rank = MPI.COMM_WORLD.rank
	f = h5py.File('stack_zp.hdf5', 'w', driver='mpio', comm=MPI.COMM_WORLD)
	print(rank)
	dset = f.create_dataset('zp_stack_100_planes', (40000, 40000, 100), chunks=(100, 100, 100), dtype='float32')

	if rank < 25:
		x = np.zeros((40000, 40000),dtype='float32')
	if rank > 74:
                x = np.zeros((40000, 40000),dtype='float32')
        else:
		x = np.array(make_zp_from_rings(250, 40000),dtype='float32')
	
	for s in range(0, 5):
		print(rank, s)
		dset.write_direct(x, source_sel=np.s_[:, :], dest_sel=np.s_[:, :, rank * 1 + s])
	f.close()
