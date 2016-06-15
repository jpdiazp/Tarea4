import scipy as sy
import matplotlib.pyplot as plt
from scipy.special import gamma
from scipy.special import gammaincinv


def chi2(x,dx):
	return (1 / (2*gamma(dx/2)) * (x/2)**(dx/2-1) * sy.exp(-x/2))
 
rand=sy.random.random(size=10000)
 
dev=120
lista1,lista2=[],[]
for i in rand:
	lista1.append(2*gammaincinv(dev/2,i))

#plt.figure("1")



x=sy.linspace(min(lista1),max(lista1),100)

for i in range(len(x)):
	lista2.append(10000*(chi2(x[i],dev)))

plt.hist(lista1, 140, facecolor='green', alpha=0.5)
plt.plot(x,lista2,'-r')	

plt.xlabel("X")
plt.ylabel("$\phi_{\mu,\sigma^2}(X)$")
plt.show()