from astropy.io import ascii
import scipy as sy
import matplotlib.pyplot as plt

data = ascii.read("datos.dat")  


#print( data["col1"])

plt.figure("1") 

plt.errorbar(data["col1"], data["col2"],data["col3"],ls="none",marker=".", color='black')
          
          
plt.xlabel("tiempo")
plt.ylabel("flujo")

#n, bins, patches = plt.hist(e, 20, facecolor='green', alpha=0.5)

x=data["col1"]
y=data["col2"]
z=data["col3"]

sig=30*(10**(-6))

Y=sy.subtract(y,1)

A=[]

for m in x:
    if m <0.4 or m>0.7:
         A.append([1,m,m**2,m**3,m**4,m**5,0]) 
  
    else:
        A.append([1,m,m**2,m**3,m**4,m**5,-1])   


theta= sy.dot( sy.linalg.inv(sy.dot( sy.transpose(A),A )) , sy.dot(sy.transpose(A),Y) )

print(theta)


plt.figure("2")

plt.errorbar(x, y,z,ls="none",marker=".", color='black')


modelo=[]

for i in x:
    e=sy.random.normal(0,sig)
    if i <0.4 or i>0.7:
        modelo.append(1+theta[0] +theta[1]*i +theta[2]*(i**2) + theta[3]*(i**3) + theta[4]*(i**4)+ theta[5]*(i**5))
  
    else:
        modelo.append(1+theta[6]*(-1)+theta[0] +theta[1]*i +theta[2]*(i**2) + theta[3]*(i**3) + theta[4]*(i**4)+ theta[5]*(i**5))   
    
    

plt.plot(x, modelo)


plt.xlabel("tiempo")
plt.ylabel("flujo")
#print(e)

plt.show()