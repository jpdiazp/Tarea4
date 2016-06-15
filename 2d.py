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
    Y=sy.subtract(y,1)    
    
    A=[]
    
    for m in x:
        if m <0.4 or m>0.7:
            A.append([1,m,m**2,m**3,m**4,m**5,0]) 
      
        else:
            A.append([1,m,m**2,m**3,m**4,m**5,-1])   
    
    theta= sy.dot( sy.linalg.inv(sy.dot( sy.transpose(A),A )) , sy.dot(sy.transpose(A),Y) )
    return theta

theta=[-0.00113006 , 0.01220323, -0.04505635 , 0.07422477 ,-0.05635102,  0.01623969,
  0.000101  ]


grados=len(x)-len(theta)

datos=[]

for l in range(1000):
    real=[]

    for i in x:
        e=sy.random.normal(0,sig)
        if i <0.4 or i>0.7:
            real.append(1+theta[0] +theta[1]*i +theta[2]*(i**2) + theta[3]*(i**3) + theta[4]*(i**4)+ theta[5]*(i**5) + e)
  
        else:
            real.append(1-theta[6]+theta[0] +theta[1]*i +theta[2]*(i**2) + theta[3]*(i**3) + theta[4]*(i**4)+ theta[5]*(i**5) + e)   
    
    
    pa=param(real)
    modelo=[]
    for s in x:
        
        if s <0.4 or s>0.7:
            modelo.append(1+pa[0] +pa[1]*s +pa[2]*(s**2) + pa[3]*(s**3) + pa[4]*(s**4)+ pa[5]*(s**5) )
  
        else:
            modelo.append(1-pa[6] +pa[0] +pa[1]*s +pa[2]*(s**2) + pa[3]*(s**3) + pa[4]*(s**4)+ pa[5]*(s**5) )   
    
    #plt.figure("3")
    #plt.plot(x,real,'ro')
    #plt.plot(x,modelo,'-b')
    #plt.show()	            
        
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
