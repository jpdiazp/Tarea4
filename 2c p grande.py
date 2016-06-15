from astropy.io import ascii
import scipy as sy
import matplotlib.pyplot as plt
from scipy.special import gammainc

data = ascii.read("datos.dat")  

x=data["col1"]
y=data["col2"]
z=data["col3"]

sig=30*(10**(-6))

def param(y):
    #Y=sy.subtract(y,1)
    Y=y    
    
    A=[]
    
    for m in x:
        A.append([1,m,m**2,m**3,m**4,m**5,m**6,m**7,m**8,m**9,m**10,m**11,m**12,m**13,m**14,m**15]) 
    
    
    theta= sy.dot( sy.linalg.inv(sy.dot( sy.transpose(A),A )) , sy.dot(sy.transpose(A),Y) )
    return theta

theta=[-0.00113006 , 0.01220323 ,-0.04505635 , 0.07422477, -0.05635102 , 0.01623969 ]


grados=len(x)-len(theta)

datos=[]

for l in range(1000):
    real=[]

    for i in x:
        e=sy.random.normal(0,sig)
        real.append(theta[0] +theta[1]*i +theta[2]*(i**2) + theta[3]*(i**3) +theta[4]*(i**4)+ theta[5]*(i**5) +e)
            
    
    pa=param(real)
    modelo=[]
    for i in x:
        modelo.append(pa[0] +pa[1]*i +pa[2]*(i**2) + pa[3]*(i**3) +pa[4]*(i**4)+pa[5]*(i**5)+pa[6]*(i**6)+pa[7]*(i**7)+pa[8]*(i**8)+pa[9]*(i**9)+pa[10]*(i**10)+pa[11]*(i**11)+pa[12]*(i**12)+pa[13]*(i**13)+pa[14]*(i**14)+pa[15]*(i**15))
        
    chi2=0
    
    for h in range(len(x)):
        
        chi2+= ((real[h]-modelo[h]) / (sig))**2
    
    #print(chi2)    
    
    pvalue=1-gammainc(grados/2,chi2/2)
    
    datos.append(pvalue)
    #print(pvalue)

plt.figure("2")

n, bins, patches = plt.hist(datos, 20, facecolor='green', alpha=0.5)


plt.xlabel("p-value")
plt.ylabel("Cantidad")


plt.show()