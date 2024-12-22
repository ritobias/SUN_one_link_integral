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

    a[0][0]=ctype(4.04623526885290664,0.74163805137086779);
    a[0][1]=ctype(0.0163846177490989801,-0.0021553588516287427);
    a[0][2]=ctype(-0.474405688605062513,-0.432057733743065147);
    a[0][3]=ctype(0.779995403479426294,-0.562037566454306991);
    a[0][4]=ctype(0.395447387176750737,0.585508785517604438);
    a[0][5]=ctype(-0.019834481701767402,0.632375717642855320);
    a[0][6]=ctype(0.443289286702868442,-0.176980401004780405);
    a[1][0]=ctype(-0.266126719511575437,0.675617015056025314);
    a[1][1]=ctype(3.94388223061273923,0.55152052487031150);
    a[1][2]=ctype(-0.487389947724269600,-0.376456000918867691);
    a[1][3]=ctype(0.303449895648367255,0.010714725439424668);
    a[1][4]=ctype(0.210113421697825853,0.402192874697527825);
    a[1][5]=ctype(-0.508151304310171583,-0.036453105528519454);
    a[1][6]=ctype(-1.055272795857492506,0.388896597392723997);
    a[2][0]=ctype(0.460112893389854394,-0.277331758517181888);
    a[2][1]=ctype(0.561421072589299377,0.057053633211441692);
    a[2][2]=ctype(4.15888892945738836,-0.06413900542213937);
    a[2][3]=ctype(-0.030833867693539181,-0.243510430022579580);
    a[2][4]=ctype(0.160472775580106874,0.159180043628188672);
    a[2][5]=ctype(0.881860172975306821,-0.193837023438362189);
    a[2][6]=ctype(0.451722993748123238,0.533536977294454829);
    a[3][0]=ctype(0.342283366745477906,0.165674106019670072);
    a[3][1]=ctype(-0.761283351672401810,0.211523506064110010);
    a[3][2]=ctype(-0.177384279723332754,0.072692415565051431);
    a[3][3]=ctype(3.80613700994691286,-0.89907312661937781);
    a[3][4]=ctype(0.77209699761009849,-1.41167195290487226);
    a[3][5]=ctype(0.755585379463069987,-0.869896507055922941);
    a[3][6]=ctype(-0.327737530140552045,-0.555741375238667217);
    a[4][0]=ctype(-0.215589406681346498,-0.211795853267406455);
    a[4][1]=ctype(0.062092124387372053,0.429675051168382581);
    a[4][2]=ctype(0.147931262316987783,-0.115414434671467227);
    a[4][3]=ctype(-1.052300189539048722,-0.887178547853267506);
    a[4][4]=ctype(4.23140638580201415,-0.82416493593448560);
    a[4][5]=ctype(-0.073275377146670544,0.558974821482671561);
    a[4][6]=ctype(0.051803631703584526,0.581605583956556034);
    a[5][0]=ctype(-0.072970863409739419,0.521025698622662581);
    a[5][1]=ctype(0.0446145374174291555,0.0228088861585820142);
    a[5][2]=ctype(-0.646324555037695866,-0.303314472483343278);
    a[5][3]=ctype(-0.458107860696342608,-0.562837439749130059);
    a[5][4]=ctype(-0.739584576279944033,0.480666030727369481);
    a[5][5]=ctype(4.25031800857740483,-0.40402892700009585);
    a[5][6]=ctype(-0.400537597449084244,0.515042929767521675);
    a[6][0]=ctype(-1.203191058106988555,-0.592658681732699319);
    a[6][1]=ctype(0.498496072497695265,0.912467541317009223);
    a[6][2]=ctype(-0.905236609028955165,0.085538508950522303);
    a[6][3]=ctype(0.593011175239783946,-0.833373642880389960);
    a[6][4]=ctype(-0.187440338623629728,0.879836412548149150);
    a[6][5]=ctype(0.442334902294961928,0.450782421485728420);
    a[6][6]=ctype(4.25812143218483623,0.48970075178061274);



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
    ftype tb=20.;
    int td=4;

    
    sun_integrator sun_int(n,td); // initialize sun_integrator object
    sun_int.set_beta(tb); // set the beta value

    ftype fsum=sun_int(a); // compute the log of the one-link integral
    //std::cout<<"fsum: "<<fsum<<", zsum: "<<std::exp(fsum)<<std::endl;
    printf("fsum: %0.9e, zsum: %0.9e\n",fsum,std::exp(fsum));

    ftype fsumr=fsum-tb*2*(td-1); // log of one-link integral after subtracting consatnt action piece
    //std::cout<<"fsumr: "<<fsumr<<", zsumr: "<<std::exp(fsumr)<<std::endl;
    printf("fsumr: %0.9e, zsumr: %0.9e\n",fsumr,std::exp(fsumr));

}