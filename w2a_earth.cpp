//
//  main.cpp
//  assgn 2
//
//  Created by 曹圳杰 on 16/9/14.
//  Copyright © 2016年 曹圳杰. All rights reserved.
//

#include <iostream>
#include <iomanip>
#include <cassert>

using namespace std;

int main()
{
    double mol = 6.02e23;
    double Earth = 5.9722e27;
    
    float Fe = 56;
    double NumberofFe = (Earth / Fe)* 26 * mol;
    double Terobyte_of_Fe = NumberofFe /8e12;
    
    float O = 16;
    double NumberofO = (Earth / O)* 8 * mol;
    double Terobyte_of_O = NumberofO /8e12;
    
    float He = 4;
    double NumberofHe = (Earth / He)* 2 * mol;
    double Terobyte_of_He = NumberofHe /8e12;
    
    cout << Terobyte_of_Fe << endl;
    cout << Terobyte_of_O << endl;
    cout << Terobyte_of_He  << endl;
    
    

    
}
