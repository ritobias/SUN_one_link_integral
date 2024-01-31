// SUN_one_link_integral.cpp : Diese Datei enthält die Funktion "main". Hier beginnt und endet die Ausführung des Programms.
//

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


    a[0][0]=ctype(5.99930869352137668,-0.01242291973879728);
    a[0][1]=ctype(-0.00612336718348499372,0.00025444820282032112);
    a[0][2]=ctype(-0.0070902446838457819,-0.0135113241418641590);
    a[0][3]=ctype(0.00181949857078487439,-0.00694246350467725905);
    a[0][4]=ctype(-0.00140879604424804475,-0.01146838389384416400);
    a[1][0]=ctype(0.00621772850776566043,0.00021252004336365737);
    a[1][1]=ctype(5.99915776058857496,0.01025247350815698);
    a[1][2]=ctype(-0.0232921099562470864,0.0174109527343487674);
    a[1][3]=ctype(0.00090727227185575634,0.01105623743148166227);
    a[1][4]=ctype(0.00075223421767953408,-0.00655854221952680970);
    a[2][0]=ctype(0.0067493572104473908,-0.0134629587255981578);
    a[2][1]=ctype(0.0231364622933113519,0.0176750152361787401);
    a[2][2]=ctype(5.99933859241426664,-0.01095968828485621);
    a[2][3]=ctype(0.00448588076610660512,-0.00555568115347733448);
    a[2][4]=ctype(-0.0178294367365768376,0.0001884045578178030);
    a[3][0]=ctype(-0.00185283642655645124,-0.00694187605225690277);
    a[3][1]=ctype(-0.00089566048174391893,0.01108660376710140597);
    a[3][2]=ctype(-0.00462503028960850397,-0.00569337002669653271);
    a[3][3]=ctype(5.99921267865563498,0.01614333964428083);
    a[3][4]=ctype(0.00092233985864315094,-0.00398962948116143719);
    a[4][0]=ctype(0.00167045094080066955,-0.01169202295946352176);
    a[4][1]=ctype(-0.00064141592666546883,-0.00662500686694504416);
    a[4][2]=ctype(0.0179042247997892880,0.0000394613868661474);
    a[4][3]=ctype(-0.00071865506352488168,-0.00402220429035570869);
    a[4][4]=ctype(5.99923240911445079,-0.00301328168921165);


 
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
    int td=4;

    sun_integrator sun_int(n,td);
    sun_int.set_beta(tb);

    ftype fsum=sun_int(a);
    std::cout<<"fsum: "<<fsum<<", zsum: "<<std::exp(fsum)<<std::endl;

    ftype fsumr=fsum-tb*2*(td-1);
    std::cout<<"fsumr: "<<fsumr<<", zsum: "<<std::exp(fsumr)<<std::endl;

}