from astropy.io import ascii
import scipy as sy
import matplotlib.pyplot as plt


data = ascii.read("datos.dat")  

x=data["col1"]
y=data["col2"]
z=data["col3"]


def mls(X,x,y,z,p):
    
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
    
    for i in X:
        
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
        

p=8

lista=[]
for i in range(1,p+1):
    lista.append(i)


k=10

Y=sy.array(y)

Y=sy.arange(Y.shape[0])

test=sy.array(sy.split(Y,len(x)/k))


medias=[]

for s in lista:
    
    
    traierror=[]    
    
    for i in test:
        trainingy=sy.array(y)
        trainingx=sy.array(x)
        
        trainingy=sy.delete(trainingy,i)
        trainingx=sy.delete(trainingx,i)
        
        model=mls(x,trainingx,trainingy,z,s)  

        error=0
        for l in i:
            error+=(Y[l]-model[0][l])**2
        
        traierror.append(error/k)

    
    media=sy.mean(sy.power(traierror,2))
    medias.append(media)
   

plt.figure("2")

plt.plot(lista, medias)

plt.xlabel("Polinomio")
plt.ylabel("CV")

plt.show()