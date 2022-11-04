def f_getRRtimeserie(mu,N):
   
    gamma=2.1
    b=0.75
    import numpy as np
    x=np.empty((N,1))
    x[:]= np.NaN
 
    k=np.empty((N,1))
    k[:]=np.NaN
 
    y=np.empty((N,1))
    y[:]=np.NaN
 
    t=np.empty((N,1))
    t[:]=np.NaN
    global tau
    tau=np.empty((N,1))
    tau[:]=np.NaN
 
    for i in range(0,N):
       
        success=False
        while not success:
            a=np.random.exponential(1)
            import random
            z=np.sign(random.random()-0.5)
            x[i,0]=np.multiply(a,z)
        
            k[i,:]=round(random.random()**(1/(-gamma+1)))
            if i<= k[i,0]:
               y[i,0]=0
            else:
             startID=int(i-k[i,0])
             endID=i
             avrg_y=np.nanmean(np.power(y[startID:endID,0],2),axis=0)
             import math
             y[i,0]= x[i,0]*math.sqrt(1+b*avrg_y)
       
            nbeats=4
            if i< nbeats+1:
               t[i,0]=0
            else:
             startID= i -np.mod(i,nbeats)-1
             endID=i-1
             t[i,0]= sum(tau[startID:endID,0])
       
            yj_sum=0
            for j in range (1,i+1):
              if j>i-k[j,0] and k[j,0]< i:
                 yj_sum=yj_sum+y[i,0]
         
            import math  
            tau[i,0] =  mu+0.03*math.sin((2*math.pi*t[i,0])/3.6)+0.025*yj_sum
        
            if tau[i,0]>2*mu or tau[i,0]<0.5*mu:
               success=False
            else:
              success=True
    #print(tau)
                
                 
                 
         
f_getRRtimeserie(752*0.001,15)           
                               
                 