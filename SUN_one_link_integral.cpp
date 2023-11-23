// SUN_one_link_integral.cpp : Diese Datei enthält die Funktion "main". Hier beginnt und endet die Ausführung des Programms.
//

#include <iostream>
#include "sun_one_link_integral.h"

int main()
{
    /*
    int n=-20;
    int nn=41;
    ftype tx=10000.5;
    bessel bessi;
    ftype bsi0;
    ftype* bsil=new ftype[nn];
    bessi(n,nn,bsil,bsi0,tx);
    for(int i=0; i<nn; ++i) {
        std::cout<<"rI("<<n+i<<","<<tx<<") = "<<bsil[i]<<std::endl;
    }
    */

    int n=5;
    ctype** a=new ctype*[n];
    for(int i=0; i<n; ++i) {
        a[i]=new ctype[n]();
        a[i][i]=0.0;
    }

    a[0][0]=ctype(5.03515381452392503,0.15076697536697338);
    a[0][1]=ctype(1.007551474092719998,0.180030690950071580);
    a[1][0]=ctype(-1.007551474092719998,0.180030690950071580);
    a[1][1]=ctype(5.03515381452392503,-0.15076697536697338);



 
    charpoly chp(n);

    ctype* cp=new ctype[n+1];
    ctype** ai=new ctype*[n];
    for(int i=0; i<n; ++i) {
        ai[i]=new ctype[n];
    }

    if(0) {
        chp(cp,ai,a);


        for(int i=0; i<=n; ++i) {
            std::cout<<cp[i]<<std::endl;
        }

        ctype res;
        for(int i=0; i<n; ++i) {
            for(int j=0; j<n; ++j) {
                res=0;
                for(int k=0; k<n; ++k) {
                    res+=a[i][k]*ai[k][j];
                }
                if(std::abs(res)>100*_fprec) {
                    std::cout<<res<<" ";
                } else {
                    std::cout<<0<<" ";
                }
            }
            std::cout<<std::endl;
        }
    }

    ftype tb=18.;

    sun_integrator sun_int(n);
    sun_int.set_beta(tb);

    ftype fsum=sun_int(a);
    std::cout<<"fsum: "<<fsum<<", zsum: "<<std::exp(fsum)<<std::endl;

}