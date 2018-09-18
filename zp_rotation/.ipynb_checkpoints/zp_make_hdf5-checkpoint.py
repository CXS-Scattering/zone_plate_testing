from mpi4py import MPI
import h5py
import numpy as np
import os

def make_zp_hdf5_file(num_slices=50,padding=50):
    rank = MPI.COMM_WORLD.rank  # The process ID (integer 0-3 for 4-process run)
    
    f = h5py.File('parallel_test.hdf5', 'w', driver='mpio', comm=MPI.COMM_WORLD)
    print(rank)
    
    os.chdir('/home/sajid/packages/zone_plate_testing/zp_rotation')
    dset = f.create_dataset('test', (40000,40000,100),chunks=(100,100,100), dtype='float32')
    
    for s in range(0,25):
        
        if rank==0 or rank==3 :
            x = np.zeros((40000,40000))
            print(rank,s,'zeros')
            dset.write_direct( x, source_sel = np.s_[:,:], dest_sel = np.s_[:,:,rank*25+s])
        else :
            os.chdir('/home/sajid/packages/zone_plate_testing/rings')
            zp = make_zp_from_rings(250,int(grid_size))
            os.chdir('/home/sajid/packages/zone_plate_testing/zp_rotation')
            print(rank,s,'zeros')
            dset.write_direct( zp, source_sel = np.s_[:,:], dest_sel = np.s_[:,:,rank*25+s])
    f.close()
    
if __name__ == "__main__":
    make_zp_hdf5_file()