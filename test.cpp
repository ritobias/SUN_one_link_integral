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

    a[0][0]=ctype(4.03153279536629126,-1.36030477486710719);
    a[0][1]=ctype(0.207635789709498481,0.852147849357878251);
    a[0][2]=ctype(-0.079167456162793984,-0.680143114079129759);
    a[0][3]=ctype(0.153810453188367999,0.684147827763677752);
    a[0][4]=ctype(-0.246373974789190770,0.191237729946770862);
    a[0][5]=ctype(-0.713703100828548864,-0.937451452749136605);
    a[0][6]=ctype(-0.0634387148746363080,0.1121534561186144406);
    a[1][0]=ctype(-0.616227595282528874,0.587122181012354661);
    a[1][1]=ctype(3.96146468044622864,1.18107156067359265);
    a[1][2]=ctype(0.401817567552687461,-0.435402494291115022);
    a[1][3]=ctype(-0.225402211409022981,0.408756082629853540);
    a[1][4]=ctype(0.267881690984102550,0.386790052358233754);
    a[1][5]=ctype(0.067751884123453644,0.408179365040684653);
    a[1][6]=ctype(-0.050256078462777712,-0.496723364042389902);
    a[2][0]=ctype(-0.704920740650506828,-0.492724848285613231);
    a[2][1]=ctype(0.183767062206161666,-0.685038700155950752);
    a[2][2]=ctype(4.07506293090021062,0.65678758790071592);
    a[2][3]=ctype(0.351677478730766432,0.004293502108527209);
    a[2][4]=ctype(0.367140568273982855,0.082086286562307437);
    a[2][5]=ctype(0.525273821930927505,-0.639637110408276924);
    a[2][6]=ctype(1.22838087078358938,1.11167469632441305);
    a[3][0]=ctype(0.020967805401555636,0.501256376801432785);
    a[3][1]=ctype(0.696269017317191510,0.560124194912971310);
    a[3][2]=ctype(-0.390458961874180407,-0.484984232799238460);
    a[3][3]=ctype(4.13295698748422681,-0.14080133597872813);
    a[3][4]=ctype(0.040848353853785647,0.394886745349070468);
    a[3][5]=ctype(-0.407178750349792361,0.540518696274436684);
    a[3][6]=ctype(-0.190013282585930588,0.014921985053688563);
    a[4][0]=ctype(0.255058113126289171,0.014054582987313786);
    a[4][1]=ctype(-0.0648608154510266296,0.0489769700608160730);
    a[4][2]=ctype(-0.418594527661728797,0.094018636177734649);
    a[4][3]=ctype(-0.325511406931579800,0.036071676069850364);
    a[4][4]=ctype(4.16082544554466354,0.92860666102667275);
    a[4][5]=ctype(0.490857300014733696,-0.345191654479777371);
    a[4][6]=ctype(0.035632027977585500,-0.271788030517615388);
    a[5][0]=ctype(-0.095400460529974518,-0.659391066599257602);
    a[5][1]=ctype(0.123118595730192436,0.321078778659311399);
    a[5][2]=ctype(-0.290873393483793916,-0.473048120547233982);
    a[5][3]=ctype(0.169713515948390847,0.263459286776934763);
    a[5][4]=ctype(-0.640465943944908828,-0.203000297867739783);
    a[5][5]=ctype(4.11647837820011661,-1.08168103643683287);
    a[5][6]=ctype(-0.289860988246427315,0.369699576869886044);
    a[6][0]=ctype(0.396781806660330482,-0.465428928767736728);
    a[6][1]=ctype(-0.322198070740070540,-0.142485726717643969);
    a[6][2]=ctype(-1.144044060704056524,0.093723865393939742);
    a[6][3]=ctype(0.240490887267577455,0.043630369472928994);
    a[6][4]=ctype(-0.445364882781047405,0.155740206315675743);
    a[6][5]=ctype(0.657695240463288708,0.345447181856773221);
    a[6][6]=ctype(4.05018298477441665,-0.48619289357478236);



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
    ftype tb=60.;
    int td=4;

    
    sun_integrator sun_int(n,td); // initialize sun_integrator object
    sun_int.set_beta(tb); // set the beta value

    ftype fsum=sun_int(a); // compute the log of the one-link integral
    //std::cout<<"fsum: "<<fsum<<", zsum: "<<std::exp(fsum)<<std::endl;
    printf("fsum: %0.8e, zsum: %0.8e\n",fsum,std::exp(fsum));

    ftype fsumr=fsum-tb*2*(td-1); // log of one-link integral after subtracting consatnt action piece
    //std::cout<<"fsumr: "<<fsumr<<", zsumr: "<<std::exp(fsumr)<<std::endl;
    printf("fsumr: %0.8e, zsumr: %0.8e\n",fsumr,std::exp(fsumr));

}