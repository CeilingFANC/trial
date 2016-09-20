# -*- coding: utf-8 -*-
"""
Created on Mon Sep 19 00:39:22 2016

@author: caozhenjie
"""

# w2c_addinghalf.py

from math import inf

s = input()


def number_from_half(s : str):
    
    num = int(s,16)
    mask_sign_16 = 1 << 15
    mask_exponent_16 = (s << 5) & 0xfffffff
    mask_fraction_16 = (s << 10) & 0xffffffff
    
    num_sign = (num & mask_sign_16)>>15
    num_exponent = (num & mask_exponent_16)& (s >> 5)
    num_fraction = (num & mask_fraction_16)& (s>> 10) 
     
    while (num_sign == 0):
        
        if num_exponent == 0:
             s = (-1)**num_sign*2**(num_exponent - 14)*(num_fraction/1024)
        elif num_exponent == 1: 
             s = inf
        else :
             s = 1 + num_fraction * (0b10)**(-10) 
        
        
    while (num_sign == 1):
        
        if num_exponent == 0:
             s = (-1)*(-1)**num_sign*2**(num_exponent - 14)*(num_fraction/1024)
        elif num_exponent == 1: 
             s = -inf
        else :
             s = (-1)*(1 + num_fraction * (0b10)**(-10)) 

    return s
    
     
def main():
    sum = 0
    sum = sum + s
    
    while True:
        try:
            sum
            
        except:
            break
    print (sum)
        


if __name__ == '__main__':
    main()
