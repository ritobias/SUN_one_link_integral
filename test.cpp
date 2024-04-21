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


    a[0][0]=ctype(2.08485652816246992,0.00024126854300312);
    a[0][1]=ctype(-0.408202909354006918,-0.165221200086766189);
    a[0][2]=ctype(-0.278898089451886340,-0.559607995600499905);
    a[0][3]=ctype(0.482071474626298121,-0.940910935653745119);
    a[0][4]=ctype(-0.343617629205371098,0.512118953356456159);
    a[0][5]=ctype(0.002208209778835291,0.969312527626935764);
    a[0][6]=ctype(-0.449423360468569406,0.403919788844845341);
    a[1][0]=ctype(-0.645340741897462149,-0.401685622573441083);
    a[1][1]=ctype(2.20921950704739287,0.53312258008536792);
    a[1][2]=ctype(-0.172555605832156863,-0.227055861708359779);
    a[1][3]=ctype(0.555279157557624138,-0.209434455829199050);
    a[1][4]=ctype(-0.250339403881251096,-0.147055682015556214);
    a[1][5]=ctype(-0.364701668964330777,0.211674594061336414);
    a[1][6]=ctype(0.445531382041857414,-0.132782761893125030);
    a[2][0]=ctype(0.405061349228038522,0.003834574582164005);
    a[2][1]=ctype(1.304825729599480467,-0.095847286033425423);
    a[2][2]=ctype(1.98106721544342949,-0.27464662138571485);
    a[2][3]=ctype(0.575129265893135461,1.188919363121698601);
    a[2][4]=ctype(0.126355134000278863,0.394386793703122433);
    a[2][5]=ctype(-0.058179074744866712,0.557271072984545550);
    a[2][6]=ctype(-0.044632094221108961,0.741534760838093896);
    a[3][0]=ctype(-0.794251947976130593,0.448932154411233280);
    a[3][1]=ctype(-0.263310314396489012,-0.328154973660015671);
    a[3][2]=ctype(-0.897862315811795761,0.183317684947314226);
    a[3][3]=ctype(1.60137556997936331,1.26515247317259940);
    a[3][4]=ctype(-0.852211436818884017,0.199758928276568188);
    a[3][5]=ctype(0.196392421504896647,-0.018511323682233974);
    a[3][6]=ctype(0.680292354842541616,-0.997496239615699697);
    a[4][0]=ctype(0.033350614210369066,0.993665944842017291);
    a[4][1]=ctype(-0.196621141220900244,0.408678006360013162);
    a[4][2]=ctype(-1.018234371377954775,-0.454615526705515045);
    a[4][3]=ctype(-0.348447050157095794,-0.135899665452070970);
    a[4][4]=ctype(1.64653208241075846,0.00505200079642452);
    a[4][5]=ctype(0.274003337497403514,-0.689781556688909968);
    a[4][6]=ctype(0.345926983260557050,-0.000044162376155817);
    a[5][0]=ctype(0.197237326495028170,0.769630864391456876);
    a[5][1]=ctype(1.205626052593284593,-0.653161789902966136);
    a[5][2]=ctype(-0.113699167849016601,-1.374502973321501367);
    a[5][3]=ctype(-0.030712545280592204,0.600738343621101871);
    a[5][4]=ctype(-0.144237391594712699,0.756086891087146396);
    a[5][5]=ctype(1.136371047476610714,0.349040078290990464);
    a[5][6]=ctype(-0.160110852233496688,0.301958340286816467);
    a[6][0]=ctype(-0.271200955634490279,0.060032004343226832);
    a[6][1]=ctype(0.542946861895871378,0.198793537618637311);
    a[6][2]=ctype(-0.087433016084948640,0.932795499765332690);
    a[6][3]=ctype(-1.52067748288727091,0.48540725554449917);
    a[6][4]=ctype(-0.144551506334741973,-0.123555050005744035);
    a[6][5]=ctype(-0.692729450318571064,-0.602292508297615877);
    a[6][6]=ctype(2.19943910228687870,0.19412675851015538);




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