import numpy as np
import subprocess 
def single_api(alpha):
	subprocess.call(["./single " + str(alpha)], shell=True)
	obj = np.loadtxt('fort.11')[0]
	dobj = np.loadtxt('fort.11')[1]
	return [obj, dobj]
	

if __name__ == "__main__":
	alpha  =1.0
	dalpha = 1e-6
	[obj1, dobj1] = single_api(alpha)
	[obj2, dobj2] = single_api(alpha + dalpha)
	fd = (obj2 - obj1)/dalpha
	print 'Finite Difference relative error', (dobj1 - fd)/obj1