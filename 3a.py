from astropy.io import ascii
import scipy as sy
import matplotlib.pyplot as plt

def mls(p):

    data = ascii.read("datos.dat")  
    
    x=data["col1"]
    y=data["col2"]
    z=data["col3"]
    
    sig=30*(10**(-6))
    
    Y=sy.subtract(y,1)
    
    A=[]

    

    for m in x:
        
        pol=[]
        for i in range(p+1):
            pol.append(m**i)
        if m <0.4 or m>0.7:
            pol.append(0)
            A.append(pol) 
      
        else:
            pol.append(-1)
            A.append(pol)
            
    
    
    theta= sy.dot( sy.linalg.inv(sy.dot( sy.transpose(A),A )) , sy.dot(sy.transpose(A),Y) )

    modelo=[]
    
    for i in x:
        
        poli=1
        for s in range(p+1):
            poli+=(theta[s]*(i**s))    
        
        e=sy.random.normal(0,sig)
        if i <0.4 or i>0.7:
            modelo.append(poli)
      
        else:
            modelo.append(poli - theta[len(theta)-1])   
        
    
    chi2=0

    for h in range(len(x)):
        
        chi2+= ((y[h]-modelo[h]) / (sig) ) **2       
        
    return modelo, theta , len(x) ,sig ,chi2
        

lista1=[]
lista2=[]

p=8

lista=[]
for i in range(1,p+1):
    lista.append(i)
    

for i in lista:
            
    datos=mls(i)

    loglike=(-datos[2]*0.5*sy.log(2*sy.pi*datos[3])) - (0.5*datos[4])
    AICc=2*len(datos[1]) - 2*loglike - ((2*len(datos[1])*(len(datos[1])-1))/(datos[2]-len(datos[1])-1))
    
    BIC=(-2)*loglike + len(datos[1])*sy.log(datos[2])
    
    lista1.append(AICc)
    lista2.append(BIC)

plt.figure("2")

print(sy.argmin(lista1) +1)
print(sy.argmin(lista2) +1)

plt.plot(lista, lista1,label="AICc")

plt.plot(lista, lista2,label="BIC")

plt.xlabel("Polinomio")
plt.ylabel("Valor")
plt.legend()

plt.show()





