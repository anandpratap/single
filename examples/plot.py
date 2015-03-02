from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)
from pylab import *

a = loadtxt("history.dat")
figure(1)
plot(a, 'rx-')
grid()
xlabel('Iteration')
ylabel(r'\alpha')
savefig('alpha.pdf')
savefig('alpha.png')

figure(2)
semilogy(abs(a-1.0), 'rx-')
grid()
xlabel('Iteration')
ylabel(r'abs(\alpha - \alpha_{\star})')
savefig('error.pdf')
savefig('error.png')

show()
