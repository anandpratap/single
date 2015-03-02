import sys
sys.path.append("../src/")
from api import single_api
import os
if __name__ == "__main__":
	f = open('history.dat', 'w')
	a0 = 2.0
	step = .5e-4
	p_benchmark = 58.2774602648079
	f.write(str(a0)+'\n')
	for i in range(30):
		[obj, dobj] = single_api(a0)
		da0 = -2*dobj*(obj - p_benchmark)*step 
		a0 = da0 + a0
		print "New a0 = ", a0
                dirname = 'run_'+str(i)
                os.system('mkdir ' + dirname)
                os.system('mv data.dat ' + dirname + '/.')
                os.system('mv fort.11 ' + dirname + '/.')
		f.write(str(a0)+'\n')
	f.close()
