# AUTHOR Fan Cao fanc@bu.edu
# w2c_addinghalf.py

from math import inf

def number_from_half(s : str):
    """return the number represented by s, a binary16 stored as a 4-character hex number"""
    g=int(s,16)
    g=(int(s,16))+int(0b10000000000000000)
    g=int(g)
    g=bin(g)
    g=str(g[3:])
    sig=int(g[:-15])
    exp=int(g[1:6],2)
    fra=g[6:]
    res=0
   

    if (exp>0)and(exp<31):
       
       exp=exp-15
       exp=2**(exp)
       fra=1+int(fra,2)*(0b10)**(-10)
       sig=(-1)**(sig)
       res=sig*exp*fra
        
    elif (exp==31)and(sig==0):
        res=inf
    
    elif(exp==31)and(sig==1):
        res=-inf
        
    else:
        exp=exp-15
        exp=2**(-14)
        fra=int(fra,2)*(0b10)**(-10)
        sig=(-1)**(sig)
        res=sig*exp*fra
    
    
        
        
        
    
    
    
    
    
    
    
    
    return res

def main():
   





try: 
   #try something
    lo=1
    savelist=[]
    for lo in range(0,999):
        
        raw=input()
        turn=int(raw,16)
        savelist.append(raw)
    
    
except: 
  #handle error
    res=0
    for z in savelist:
        a=number_from_half(z)
        
        
        res=res+a
    print(res)
    
print("spotify have too much advertise, and it is not that good!")

print("how can I see the changes in all those files")

print (" I would like to have some cream")

print( "I SHOULD DELETE SOMETHING")





    # hint1: use int(input(),16)
    # hint2: use try: except: to halt


if __name__ == '__main__':
    main()
