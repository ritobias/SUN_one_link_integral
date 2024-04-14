/*
    cpp file for testing the sun_integrator class
    Copyright (C) 2024  Tobias Rindlisbacher

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.

    Email: ritobias@gmx.ch
*/
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

    int n=7;
    ctype** a=new ctype*[n];
    for(int i=0; i<n; ++i) {
        a[i]=new ctype[n]();
        a[i][i]=0.0;
    }


    a[0][0]=ctype(1.77032257780899490,1.17637090010593349);
    a[0][1]=ctype(-0.267962957024662790,0.084184023899804321);
    a[0][2]=ctype(0.278381745076954550,-0.660987989095463491);
    a[0][3]=ctype(-1.145184420735959282,-0.036418247090869084);
    a[0][4]=ctype(-0.499876741602245023,0.023956729771446555);
    a[0][5]=ctype(0.823186233861974775,-0.716811348086584297);
    a[0][6]=ctype(0.074515860530701684,0.141757010130999301);
    a[1][0]=ctype(0.537636878050848206,-0.194320841176854463);
    a[1][1]=ctype(1.56586164484456644,0.98282601168938169);
    a[1][2]=ctype(0.486647001406682026,-0.100249941335183094);
    a[1][3]=ctype(0.273125354661653391,1.101560635183096207);
    a[1][4]=ctype(-0.156933715325590274,-0.895163598979136522);
    a[1][5]=ctype(-0.033151878800381148,-1.048942323553379972);
    a[1][6]=ctype(-0.962102019549421243,-0.332760044830958647);
    a[2][0]=ctype(-0.024572392394026647,-1.126699130651564018);
    a[2][1]=ctype(-0.582091015929821374,-0.774416001588208849);
    a[2][2]=ctype(1.337036671223019306,-0.094444756865582703);
    a[2][3]=ctype(0.733277499727362727,-0.319578635758831560);
    a[2][4]=ctype(-1.025308212281614987,0.049250632427100148);
    a[2][5]=ctype(0.119843897066049520,-0.404095947749221307);
    a[2][6]=ctype(1.037148188621221220,-0.180513655970973694);
    a[3][0]=ctype(-0.083047447807247077,0.820406755264060836);
    a[3][1]=ctype(-1.01587074374443446,1.19639644509523047);
    a[3][2]=ctype(0.133116531274068275,-0.125354462506132124);
    a[3][3]=ctype(1.85307967300817528,0.61773894794579203);
    a[3][4]=ctype(0.04234017460242573,-1.55099295207656419);
    a[3][5]=ctype(0.562486698912928965,-0.027818762395531628);
    a[3][6]=ctype(0.728332947353406878,-0.394698225721385149);
    a[4][0]=ctype(0.72706560545735989,-1.27109781276744006);
    a[4][1]=ctype(-0.746928613991926655,-0.905903183702451187);
    a[4][2]=ctype(0.781381688728521827,-0.489384449062277983);
    a[4][3]=ctype(-0.894293851251272203,-0.740270956550330789);
    a[4][4]=ctype(1.73459713839945853,-0.88143539433324801);
    a[4][5]=ctype(-0.562038424592954781,0.751989372628345038);
    a[4][6]=ctype(0.145453407084046489,0.717834048635149790);
    a[5][0]=ctype(-0.382076831393577241,-0.652306289739303186);
    a[5][1]=ctype(-0.586002508032405855,-1.043911696254844198);
    a[5][2]=ctype(0.753391156128081530,0.043876340251929956);
    a[5][3]=ctype(-0.113671207108066586,0.519758058487388984);
    a[5][4]=ctype(0.292489204284583949,0.455911560133112822);
    a[5][5]=ctype(1.57157436472616021,-0.19305909904857020);
    a[5][6]=ctype(-1.87658126784836152,-0.16300197418584520);
    a[6][0]=ctype(0.477358273513711243,0.139171138293549127);
    a[6][1]=ctype(-0.247389743883350174,0.464068449996036216);
    a[6][2]=ctype(-0.0542173491691491065,0.0806776833155646109);
    a[6][3]=ctype(-1.94346862000794652,1.11098972810968510);
    a[6][4]=ctype(1.14698484018323441,0.91356279553067140);
    a[6][5]=ctype(-0.034777331388264657,-0.509156754323076433);
    a[6][6]=ctype(2.12111669278155353,-0.31216750946540253);



      if(0) {
        // test charpoly:
        charpoly chp(n);
        ctype* cp=new ctype[n+1];
        ctype** ai=new ctype*[n];
        for(int i=0; i<n; ++i) {
            ai[i]=new ctype[n];
        }

        chp(cp,ai,a); // get the characteristic polynomial coefficients cp[] and the matrix inverse ai[][] of a[][]

        for(int i=0; i<=n; ++i) {
            std::cout<<cp[i]<<std::endl;
        }

        // check accuracy of inverse ai[][] of a[][]:
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

        delete[] cp;
        for(int i=0; i<n; ++i) {
            delete[] ai[i];
        }
        delete[] ai;
    }

    // specifiy a beta-value and number of dimensions:
    ftype tb=21.;
    int td=4;

    
    sun_integrator sun_int(n,td); // initialize sun_integrator object
    sun_int.set_beta(tb); // set the beta value

    ftype fsum=sun_int(a); // compute the log of the one-link integral
    std::cout<<"fsum: "<<fsum<<", zsum: "<<std::exp(fsum)<<std::endl;

    ftype fsumr=fsum-tb*2*(td-1); // log of one-link integral after subtracting consatnt action piece
    std::cout<<"fsumr: "<<fsumr<<", zsum: "<<std::exp(fsumr)<<std::endl;

}