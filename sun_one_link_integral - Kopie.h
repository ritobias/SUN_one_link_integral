#pragma once
#include <algorithm>
#include <complex>
#include <limits>
#include <cmath>

typedef double ftype;
typedef std::complex<ftype> ctype;

static const ftype _fmin=std::numeric_limits<ftype>::min();
static const ftype _fmax=std::numeric_limits<ftype>::max();
static const ftype _fprec=std::numeric_limits<ftype>::epsilon();

namespace levi_civita {
	static void permute(int** tlcs,int* arr,const int& arrlen,int& icnt,int sign=1,int ipos=0) {
		//recursive function that generates all permutations of 0 1 2 3 ... d 
		//and the corresponding signs {-1,1}.

		int tval=0;
		if(ipos==arrlen-1) {
			for(int i=0; i<arrlen; ++i) {
				tlcs[icnt][i]=arr[i];
			}
			tlcs[icnt][arrlen]=sign;
			++icnt;
		} else {
			permute(tlcs,arr,arrlen,icnt,sign,ipos+1);
			for(int i=ipos+1; i<arrlen; ++i) {
				tval=arr[ipos];
				arr[ipos]=arr[i];
				arr[i]=tval;
				permute(tlcs,arr,arrlen,icnt,-sign,ipos+1);
				arr[i]=arr[ipos];
				arr[ipos]=tval;
			}
		}
	}

	static void new_lcs(int n,int** lcs, int& dfac) {
		//sets up lcs[][], such that lcs[i][j] contains for "j<d" the "j"-th coordinate of the "i"-th non-zero entry
		//of the Levi-Civita symbol (assuming that they are ordered such that \epsilon_{0 1 2 3} is the first, and
		//and \epsilon_{3 2 1 0} the last nonzero entry) for "j=d" its value {-1,1}. 
		//For example \epsilon_{0 1 2 3} = 1 is represented by lcs[0],
		//where: lcs[0][0]=0, lcs[0][1]=1, lcs[0][2]=2, lcs[0][3]=3, lcs[0][4]=1, and
		//as another example: \epsilon_{0 1 3 2} = -1 is represented by lcs[1],
		//where: lcs[1][0]=0, lcs[1][1]=1, lcs[1][2]=3, lcs[1][3]=2, lcs[1][4]=-1. 

		dfac=1;
		for(int i=2; i<=n; ++i) {
			dfac*=i;
		}
		lcs=new int* [dfac];
		for(int i=0; i<dfac; ++i) {
			lcs[i]=new int[n+1];
		}

		int icnt=0;
		int* tarr;
		tarr=new int[n];
		for(int i=0; i<n; ++i) {
			tarr[i]=i;
		}
		permute(lcs,tarr,n,icnt);
		delete[] tarr;
		if(0) { //set to 1 for printing the entries of lcs[][]
			for(int i=0; i<dfac; ++i) {
				for(int j=0; j<n; ++j) {
					std::cout<<lcs[i][j]<<" ";
				}
				std::cout<<": "<<lcs[i][n]<<std::endl;
			}
		}

	}

	static void delete_lcs(int n,int** lcs) {
		int dfac=1;
		for(int i=2; i<=n; ++i) {
			dfac*=i;
		}
		for(int i=0; i<dfac; ++i) {
			delete[] lcs[i];
		}
		delete[] lcs;
	}
};

class charpoly {
public:
	charpoly(): n(0),lb(0)	{

	}

	charpoly(int tn): n(tn) {
		lb=new ctype**[2];
		for(int i=0; i<2; ++i) {
			new_matrix(lb[i],0);
		}
	}

	void set_n(int tn) {
		if(tn>0) {
			if(tn!=n) {
				n=tn;
				if(lb!=0) {
					for(int i=0; i<2; ++i) {
						delete_matrix(lb[i]);
					}
					delete[] lb;
				}
				lb=new ctype**[2];
				for(int i=0; i<2; ++i) {
					new_matrix(lb[i],0);
				}
			}
		} else {
			n=0;
			if(lb!=0) {
				for(int i=0; i<2; ++i) {
					delete_matrix(lb[i]);
				}
				delete[] lb;
				lb=0;
			}
		}
	}

	~charpoly() {
		for(int i=0; i<2; ++i) {
			delete_matrix(lb[i]);
		}
		delete[] lb;
	}

	void operator()(ctype* lc,ctype** ta) {
		int i,i1,i2,i3;
		int ilb=0;
		lc[n]=1.0;

		ctype tlc,ttlc;
		ctype** tlb;
		ctype** tlbo;
		ctype* ta1;
		ctype* tlb1;

		tlb=lb[ilb];
		tlc=0;
		ttlc=-lc[n];
		for(i1=0; i1<n; ++i1) {
			ta1=ta[i1];
			tlb1=tlb[i1];
			for(i2=0; i2<n; ++i2) {
				if(i1==i2) {
					tlb1[i2]=ttlc;
					tlc+=ttlc*ta1[i2];
				} else {
					tlb1[i2]=0;
				}
			}
		}
		lc[n-1]=tlc;
		ilb=1-ilb;
		

		for(i=1; i<n; ++i) {
			tlb=lb[ilb];
			tlbo=lb[1-ilb];
			tlc=0;
			ttlc=-lc[n-i];
			for(i1=0; i1<n; ++i1) {
				ta1=ta[i1];
				tlb1=tlb[i1];
				for(i2=0; i2<n; ++i2) {
					if(i1==i2) {
						tlb1[i2]=ttlc;
					} else {
						tlb1[i2]=0;
					}
					for(i3=0; i3<n; ++i3) {
						tlb1[i2]+=ta1[i3]*tlbo[i3][i2];
					}
					tlc+=tlb1[i2]*ta[i2][i1];
				}
			}
			lc[n-i-1]=tlc/(ftype)(i+1);
			ilb=1-ilb;
		}
	}

	void operator()(ftype* lc,ctype** ta) {
		// assuming real charactersitic polynomial coefficients (from hermitian ta[i][j])
		int i,i1,i2,i3;
		int ilb=0;
		lc[n]=1.0;

		ftype tlc;
		ftype ttlc;
		ctype** tlb;
		ctype** tlbo;
		ctype* ta1;
		ctype* tlb1;

		tlb=lb[ilb];
		tlc=0;
		ttlc=-lc[n];
		for(i1=0; i1<n; ++i1) {
			ta1=ta[i1];
			tlb1=tlb[i1];
			for(i2=0; i2<n; ++i2) {
				if(i1==i2) {
					tlb1[i2]=ttlc;
					tlc+=ttlc*std::real(ta1[i2]);
				} else {
					tlb1[i2]=0;
				}
			}
		}
		lc[n-1]=tlc;
		ilb=1-ilb;


		for(i=1; i<n; ++i) {
			tlb=lb[ilb];
			tlbo=lb[1-ilb];
			tlc=0;
			ttlc=-lc[n-i];
			for(i1=0; i1<n; ++i1) {
				ta1=ta[i1];
				tlb1=tlb[i1];
				for(i2=0; i2<n; ++i2) {
					if(i1==i2) {
						tlb1[i2]=ttlc;
					} else {
						tlb1[i2]=0;
					}
					for(i3=0; i3<n; ++i3) {
						tlb1[i2]+=ta1[i3]*tlbo[i3][i2];
					}
					tlc+=std::real(tlb1[i2]*std::conj(ta1[i2]));
				}
			}
			lc[n-i-1]=tlc/(ftype)(i+1);
			ilb=1-ilb;
		}
	}

private:
	void new_matrix(ctype**& ta,int init=-1) {
		ta=new ctype*[n];
		if(init==-1) {
			for(int ic1=0; ic1<n; ++ic1) {
				ta[ic1]=new ctype[n];
			}
		} else if(init==0) {
			for(int ic1=0; ic1<n; ++ic1) {
				ta[ic1]=new ctype[n]();
			}
		} else {
			for(int ic1=0; ic1<n; ++ic1) {
				ta[ic1]=new ctype[n]();
				ta[ic1][ic1]=1;
			}
		}
	}

	void delete_matrix(ctype**& ta) {
		for(int ic1=0; ic1<n; ++ic1) {
			delete[] ta[ic1];
		}
		delete[] ta;
		ta=0;
	}

	int n;
	ctype*** lb;
};

class bessel {
public:
	bessel() {
		iacc=40;
		bigno=1e10;
		bigni=1e-10;

		p01=1.0; p02=3.5156229; p03=3.0899424; p04=1.2067492;
		p05=0.2659732; p06=0.360768e-1; p07=0.45813e-2;
		q01=0.39894228; q02=0.1328592e-1; q03=0.225319e-2;
		q04=-0.157565e-2; q05=0.916281e-2; q06=-0.2057706e-1;
		q07=0.2635537e-1; q08=-0.1647633e-1; q09=0.392377e-2;

		p11=0.5; p12=0.87890594; p13=0.51498869; p14=0.15084934;
		p15=0.2658733e-1; p16=0.301532e-2; p17=0.32411e-3;
		q11=0.39894228; q12=-0.3988024e-1; q13=-0.362018e-2;
		q14=0.163801e-2; q15=-0.1031555e-1; q16=0.2282967e-1;
		q17=-0.2895312e-1; q18=0.1787654e-1; q19=-0.420059e-2;
	}

	~bessel() {

	}

	ftype operator()(int n,ftype x) {
		/*----------------------------------------------------------------------
		!     this subroutine calculates the first kind modified bessel function
		!     of integer order n, for any real x.
		!     it uses continued fractions (cf. https://dlmf.nist.gov/10.33.1)
		!     which are computed in terms of "continuants".
		------------------------------------------------------------------------*/

		if(n<0) n=-n;
		if(n==0)  return (bessi0(x));
		if(n==1)  return (bessi1(x));
		if(x==0) return 0;

		ftype tox=2.0/x;
		ftype bip=0;
		ftype bi=1.0;
		ftype bsi=0;
		ftype bim=0;
		int m=2*(n+(int)std::floor(std::sqrt(iacc*n)));
		int m2=(int)std::ceil(0.5*(1.0+std::sqrt(1.0+iacc*x*x)));
		if(m<m2) {
			m=m2;
		}
		for(int j=m; j>0; j--) {
			bim=bip+j*tox*bi;
			bip=bi;
			bi=bim;
			if(std::abs(bi)>bigno) {
				bi*=bigni;
				bip*=bigni;
				bsi*=bigni;
			}
			if(j==n)  bsi=bip;
		}
		return (bsi*bessi0(x)/bi);
	}

	void operator()(int n,int nn,ftype* bsil,ftype& bsi0,ftype x) {
		/*----------------------------------------------------------------------
		!     this subroutine calculates the first kind modified bessel function
		!     of integer orders [n,...n+nn-1], for any real x. 
		!     it uses continued fractions (cf. https://dlmf.nist.gov/10.33.1)
		!     which are computed in terms of "continuants".
		------------------------------------------------------------------------*/
		int tn0,tn1,revo;
		int ttn1=n+nn-1;
		int ttn0=n;
		if(std::abs(ttn0)>std::abs(ttn1)) {
			tn0=-ttn1;
			tn1=-ttn0;
			revo=1;
		} else {
			tn0=ttn0;
			tn1=ttn1;
			revo=0;
		}
		if(tn0<=0) {
			tn0=1;
		}

		ftype tox=2.0/x;
		ftype bip=0;
		ftype bi=1.0;
		ftype bim=0;
		ftype bsi=0;
		int i,j,k,i0;
		j=2*(tn1+(int)std::floor(std::sqrt(iacc*(ftype)tn1)));
		i=(int)std::ceil(0.5*(1.0+std::sqrt(1.0+iacc*x*x)));
		if(j<i) {
			j=i;
		}
		while(j>tn1) {
			bim=bip+j*tox*bi;
			bip=bi;
			bi=bim;
			if(std::abs(bi)>bigno) {
				bi*=bigni;
				bip*=bigni;
			}
			--j;
		}

		if(revo==0) {
			i=nn-1;
			while(j>=tn0) {
				bim=bip+j*tox*bi;
				bip=bi;
				bi=bim;
				if(std::abs(bi)>bigno) {
					bi*=bigni;
					bip*=bigni;
					for(k=i+1; k<nn; ++k) {
						bsil[k]*=bigni;
					}
				}
				bsil[i]=bip;
				--j;
				--i;
			}
			if(i>=0) {
				i0=i;
				bsil[i0]=bi;
				--i;
				while(i>=0) {
					bsil[i]=bsil[2*i0-i];
					--i;
				}
			}
		} else {
			i=0;
			while(j>=tn0) {
				bim=bip+j*tox*bi;
				bip=bi;
				bi=bim;
				if(std::abs(bi)>bigno) {
					bi*=bigni;
					bip*=bigni;
					for(k=0; k<i; ++k) {
						bsil[k]*=bigni;
					}
				}
				bsil[i]=bip;
				--j;
				++i;
			}
			if(i<nn) {
				i0=i;
				bsil[i0]=bi;
				++i;
				while(i<nn) {
					bsil[i]=bsil[2*i0-i];
					++i;
				}
			}
		}

		while(j>0) {
			bim=bip+j*tox*bi;
			bip=bi;
			bi=bim;
			if(bi>bigno) {
				bi*=bigni;
				bip*=bigni;
				for(k=0; k<nn; ++k) {
					bsil[k]*=bigni;
				}
			}
			--j;
		}
		bsi=1.0/bi;
		for(k=0; k<nn; ++k) {
			bsil[k]*=bsi;
		}
		bsi0=bessi0(x);
	}

	ftype bessi0(ftype x) {
		ftype y,ax;
		ax=std::abs(x);
		if(ax<3.75) {
			y=(x/3.75);
			y*=y;
			return (p01+y*(p02+y*(p03+y*(p04+y*(p05+y*(p06+y*p07))))));
		} else {
			ftype bx;
			y=3.75/ax;
			bx=std::exp(ax)/std::sqrt(ax);
			ax=q01+y*(q02+y*(q03+y*(q04+y*(q05+y*(q06+y*(q07+y*(q08+y*q09)))))));
			return (ax*bx);
		}
	}

	ftype bessi1(ftype x) {
		ftype y,ax;
		ax=std::abs(x);
		if(ax<3.75) {
			y=(x/3.75);
			y*=y;
			return(x*(p11+y*(p12+y*(p13+y*(p14+y*(p15+y*(p16+y*p17)))))));
		} else {
			ftype bx;
			y=3.75/ax;
			bx=std::exp(ax)/std::sqrt(ax);
			ax=q11+y*(q12+y*(q13+y*(q14+y*(q15+y*(q16+y*(q17+y*(q18+y*q19)))))));
			if(x>0) {
				return (ax*bx);
			} else {
				return -(ax*bx);
			}
		}
	}

	int iacc;
	ftype bigno;
	ftype bigni;
	ftype p01,p02,p03,p04,p05,p06,p07;
	ftype q01,q02,q03,q04,q05,q06,q07,q08,q09;
	ftype p11,p12,p13,p14,p15,p16,p17;
	ftype q11,q12,q13,q14,q15,q16,q17,q18,q19;
};

class esolve_h {
public:
	esolve_h(): n(0),d(0),a(0),s(0),w(0) {

	}

	esolve_h(int tn): n(tn),d(0),a(0) {
		s=new ftype[n-1]();
		w=new ctype[n];
	}

	void set_n(int tn) {
		if(tn>0) {
			if(tn!=n) {
				n=tn;
				d=0;
				a=0;
				if(s!=0) {
					delete[] s;
				}
				s=new ftype[n-1]();
				if(w!=0) {
					delete[] w;
				}
				w=new ctype[n];
			}
		} else {
			n=0;
			d=0;
			a=0;
			if(s!=0) {
				delete[] s;
				s=0;
			}
			if(w!=0) {
				delete[] w;
				w=0;
			}
		}
	}

	~esolve_h() {
		if(s!=0) {
			delete[] s;
		}
		if(w!=0) {
			delete[] w;
		}
	}

	void to_tridiag() {
		// Housholder reduction of hermitian nxn matrix "a" to tri-diagonal
		// form. The n real diagonal entries ares stored in the array "d" and
		// the n-1 real sub-/super-diagonal entries are stored in the array "s".
		// Upon exit {a[j][i]}_{j=i+1,...,n} is the i-th complex vector u_i
		// and the diagonal element a[i][i] is the i-th complex scalar g_i
		// so that (q_i)_{ab}=\delta_{ab} - conj(g_i) (u_i)_a (conj(u_i))_b
		// is the i-th Housholder reflector. The product of the n-1 reflectos
		// {q_i}_{i=0,...,n-1} yields the U(n) matrix Q=q_{n-2} q_{n-1} ... q_{0}
		// for which  Q.a.Q^{\dagger} = T, where "a" refers to the original
		// matrix "a" and T to the tri-diagonal matrix with diagonal "d" and
		// sub-/super-diagonal "s".

		ftype b,bi,t;
		int i,j,k;
		ctype al,e,ei,g;
		for(k=0; k<n-1; ++k) {
			b=0;
			for(i=k+1; i<n; ++i) {
				bi=std::abs(a[i][k]);
				if(bi>b) {
					b=bi;
				}
			}
			g=0;
			if(b>0) {
				bi=1/b;
				t=0;
				for(j=k+1; j<n; ++j) {
					ctype& ta=a[j][k];
					ta*=bi;
					t+=std::norm(ta);
				}
				t=std::sqrt(t);
				if(a[k+1][k].real()<0) {
					t=-t;
				}
				e=a[k+1][k]+t;
				ei=1.0/e;
				a[k+1][k]=1.0;
				for(j=k+2; j<n; ++j) {
					a[j][k]*=ei;
				}
				g=e/t;
				t*=b;

				for(j=k+1; j<n; ++j) {
					w[j]=0;
					for(i=j+1; i<n; ++i) {
						w[j]+=std::conj(a[i][j])*a[i][k];
					}
				}

				al=0;
				for(j=k+1; j<n; ++j) {
					ctype& ta=a[j][k];
					for(i=j; i<n; ++i) {
						w[i]+=a[i][j]*ta;
					}
					w[j]*=-g;
					al+=std::conj(w[j])*ta;
				}
				al*=-0.5*g;
				for(j=k+1; j<n; ++j) {
					w[j]+=al*a[j][k];
				}
				for(j=k+1; j<n; ++j) {
					for(i=j; i<n; ++i) {
						a[i][j]+=w[i]*std::conj(a[j][k])+a[i][k]*std::conj(w[j]);
					}
				}
				s[k]=-t;
			}
			d[k]=a[k][k].real();
			a[k][k]=g;
		}
		d[n-1]=a[n-1][n-1].real();
		a[n-1][n-1]=0;
	}

	void mq(ftype& e1,ftype& e2,const ftype& d1,const ftype& d2,const ftype& s1) {
		// modified quadratic formula for the eigenvalues of a real 2x2 matrix.
		// As we consider here only hermitian matrices, there should be no
		// complex eigenvalues. We therefore assume that "disc<0" occurse only
		// due to roundoff errors and should be interpreted as disc=0.

		ftype t=0.5*(d1+d2);
		ftype dt=d1*d2-s1*s1;
		ftype disc=t*t-dt;
		if(disc<0) {
			disc=0;
		} else {
			disc=std::sqrt(disc);
		}
		if(t>0) {
			if(t+disc>0) {
				e1=t+disc;
				e2=dt/e1;
			} else {
				e1=0;
				e2=0;
			}
		} else {
			if(t-disc<0) {
				e1=t-disc;
				e2=dt/e1;
			} else {
				e1=0;
				e2=0;
			}
		}
	}

	void uvtrid(ftype& v1,ftype& v2,ftype& v3,int i,int f) {
		// computes the first column of the double-shifted (f-i)x(f-i) sub-matrix.
		// (cf. Ex. 5.7.40, in "D.S. Watkins - Fundamentals of Matrix Computations",
		//  2nd Ed., ISBN 0-471-21394-2)
		v1=((d[i]-d[f-1])*(d[i]-d[f])-s[f-1]*s[f-1])/s[i]+s[i];
		v2=d[i]+d[i+1]-d[f-1]-d[f];
		if(i+1<n-1) {
			v3=s[i+1];
		} else {
			v3=0;
		}
	}

	void uqtrid(ftype& g,ftype& u2,ftype& u3,const ftype& y1,const ftype& y2, const ftype& y3) {
		// computes reflector for double-shifted QR iteration
		// (cf. Ex. 5.7.41, in "D.S. Watkins - Fundamentals of Matrix Computations",
		//  2nd Ed., ISBN 0-471-21394-2)
		ftype t=std::sqrt(y1*y1+y2*y2+y3*y3);
		if(y1<0) {
			t=-t;
		}
		if(std::abs(t)>0) {
			g=(y1+t)/t;
			u2=y2/(y1+t);
			u3=y3/(y1+t);
		} else {
			g=0;
			u2=0;
			u3=0;
		}
	}

	void initblg(ftype& b1,ftype& b2,ftype& b3,int i, const ftype& y1, const ftype& y2, const ftype& y3) {
		// initializes a QR iteration by inserting a "bulge" (b1,b2,b3) to the lower-right submatrix
		// that starts at the i-th row and column.
		// (cf. Sec. 5.7 "The Double-Step QR Algorithm"
		//  in "D.S. Watkins - Fundamentals of Matrix Computations", 2nd Ed., ISBN 0-471-21394-2)
		// Adaptet to operate on tri-diagonal matrix with diagonal "d" and sub-/super-diagonal "s".
		ftype g,u2,u3;
		uqtrid(g,u2,u3,y1,y2,y3);
		ftype gsq=g*g;
		ftype a[4];
		a[0]=d[i]+u2*s[i];
		a[1]=s[i];
		if(i+1<n) {
			a[1]+=u2*d[i+1];
			if(i+2<n) {
				a[1]+=u3*s[i+1];
				a[2]=u2*s[i+1]+u3*d[i+2];
			} else {
				a[2]=0;
			}
			if(i+3<n) {
				a[3]=u3*s[i+2];
			} else {
				a[3]=0;
			}
		} else {
			a[2]=0;
			a[3]=0;
		}

		ftype au=a[0]+a[1]*u2+a[2]*u3;

		d[i]+=gsq*au-2.0*g*a[0];
		if(i+1<n) {
			d[i+1]+=gsq*au*u2*u2-2.0*g*a[1]*u2;
			s[i]+=gsq*au*u2-g*(a[1]+u2*a[0]);

			if(i+2<n) {
				d[i+2]+=gsq*au*u3*u3-2.0*g*a[2]*u3;
				s[i+1]+=gsq*au*u3*u2-g*(a[2]*u2+u3*a[1]);

				if(i+3<n) {
					s[i+2]-=g*a[3]*u3;
				}
			}
		}
		b1=gsq*au*u3-g*(a[2]+u3*a[0]);
		b2=-g*a[3];
		b3=b2*u2;
	}

	void chaseblg(ftype& b1,ftype& b2,ftype& b3,int k) {
		// Advances the QR iteration by moving the "bulge" form the k-th to the (k+1)-th
		// column. (cf. Sec. 5.7 "The Double-Step QR Algorithm"
		//  in "D.S. Watkins - Fundamentals of Matrix Computations", 2nd Ed., ISBN 0-471-21394-2)
		// Adaptet to operate on tri-diagonal matrix with diagonal "d" and sub-/super-diagonal "s".

		ftype g,u2,u3;
		uqtrid(g,u2,u3,s[k-1],b1,b2);
		ftype gsq=g*g;
		ftype a[5];
		a[0]=s[k-1]+u2*b1+u3*b2;
		a[1]=d[k]+u2*s[k]+u3*b3;
		a[2]=s[k];
		a[3]=b3;
		if(k+1<n) {
			a[2]+=u2*d[k+1];
			if(k+2<n) {
				a[2]+=u3*s[k+1];
				a[3]+=u2*s[k+1]+u3*d[k+2];
				if(k+3<n) {
					a[4]=u3*s[k+2];
				} else {
					a[4]=0;
				}
			} else {
				a[4]=0;
			}
		} else {
			a[4]=0;
		}
		ftype au=a[1]+a[2]*u2+a[3]*u3;

		d[k]+=gsq*au-2.0*g*a[1];
		s[k-1]-=g*a[0];
		if(k+1<n) {
			d[k+1]+=gsq*au*u2*u2-2.0*g*u2*a[2];
			s[k]+=gsq*au*u2-g*(u2*a[1]+a[2]);
			if(k+2<n) {
				d[k+2]+=gsq*au*u3*u3-2.0*g*u3*a[3];
				s[k+1]+=gsq*au*u3*u2-g*(u3*a[2]+a[3]*u2);
				if(k+3<n) {
					s[k+2]-=g*a[4]*u3;
				}
			}
		}
		b1=b3+gsq*au*u3-g*(u3*a[1]+a[3]);
		b2=-g*a[4];
		b3=b2*u2;
	}

	int evs_from_tridiag() {
		// implements the implicitly shifted double-step QR algorithm for real, tri-diagonal
		// nxn matrices, given in terms of their n diagonal entries in "d" and their (n-1)
		// sub-/super-diagonal entries in "s". Upon exit the n real eigenvalues are contained
		// in the array "d" and the entries of the array "s" are set to zero.
		// (cf.Sec. 5.7 "The Double-Step QR Algorithm"
		//  in "D.S. Watkins - Fundamentals of Matrix Computations", 2nd Ed., ISBN 0-471-21394-2)

		int ek=n-1;
		int sk;
		int nit=0;
		int nitmax=10*n; // maximum number of QR iterations between the finding of eigenvalues
		                 // (10*n is probably way too large)
		int l;
		ftype y1,y2,y3,b1,b2,b3;
		while(ek>0&&nit<nitmax) {
			sk=ek;
			// search for lowest laying diagonal block (deflation):
			while(sk>0) {
				if(std::abs(s[sk-1])<_fprec*(std::abs(d[sk-1])+std::abs(d[sk]))) {
					s[sk-1]=0;
					break;
				} else {
					--sk;
				}
			}
			if(sk==ek) {
				// if 1x1 block: d[sk] contains eigenvalue;
				// continue finding next block;
				ek=sk-1;
				nit=0;
			} else if(sk==ek-1) {
				// if 2x2 block: determine eigenvalues with modified quadratic formula "mq";
				mq(y1,y2,d[sk],d[ek],s[sk]);
				// set d[sk], d[ek] to the values of the found eigenvalues and set s[sk] to zero;
				d[sk]=y1;
				d[ek]=y2;
				s[sk]=0;
				// continue finding next block;
				ek=sk-1;
				nit=0;
			} else {
				//if block is bigger than 2x2: perform QR iteration on block
				uvtrid(y1,y2,y3,sk,ek);
				initblg(b1,b2,b3,sk,y1,y2,y3);
				for(l=sk+1; l<ek; ++l) {
					chaseblg(b1,b2,b3,l);
				}
				//continue to serach for lowest laying diagonal block in updated matrx;
			}
			++nit;
		}
		std::sort(d,d+n,std::greater<ftype>());
		if(nit==nitmax) {
			return 0;
		} else {
			return 1;
		}
	}

	int operator()(ftype* td,ctype** ta) {
		a=ta;
		d=td;
		to_tridiag();
		return evs_from_tridiag();
	}

private:
	int n;
	ftype* s;
	ftype* d;
	ctype** a;
	ctype* w;
};

class sun_integrator {
public:
	sun_integrator(int tn): n(tn),es(tn),bess(),lcs(0),nlcs(0) {
		levi_civita::new_lcs(n,lcs,nlcs);
		new_matrix(a,0);
		d=new ftype[n]();
		m=new ftype[n]();
		eqprec=0.00001;
	}

	~sun_integrator() {
		levi_civita::delete_lcs(n,lcs);
		delete_matrix(a);
		delete[] d;
		delete[] m;
	}

	ftype operator()(ftype beta, ctype** staple) {
		ftype trs=beta/(ftype)n;
		trs*=trs;
		matrix_mult_an(staple,staple,trs,a);
		es(d,a);
		int i,j;
		m[0]=1;
		for(i=1; i<n; ++i) {
			if(std::abs(d[i-1]-d[i])<0.5*eqprec*(std::abs(d[i-1])+std::abs(d[i]))) {
				m[i]=m[i-1]+1;
			} else {
				m[i]=1;
			}
		}

		ftype lvdm=0;
		ftype tdi,tli;
		for(i=0; i<n-1; ++i) {
			tdi=d[i];
			tli=std::log((ftype)(1+i));
			for(j=i+1; j<n; ++j) {
				if(m[j]==1) {
					lvdm+=std::log(std::abs(tdi-d[j]))-tli;
				}
			}
		}

	}

private:
	void new_matrix(ctype**& ta,int init=-1) {
		ta=new ctype*[n];
		if(init==-1) {
			for(int ic1=0; ic1<n; ++ic1) {
				ta[ic1]=new ctype[n];
			}
		} else if(init==0) {
			for(int ic1=0; ic1<n; ++ic1) {
				ta[ic1]=new ctype[n]();
			}
		} else {
			for(int ic1=0; ic1<n; ++ic1) {
				ta[ic1]=new ctype[n]();
				ta[ic1][ic1]=1;
			}
		}
	}

	void delete_matrix(ctype**& ta) {
		for(int ic1=0; ic1<n; ++ic1) {
			delete[] ta[ic1];
		}
		delete[] ta;
		ta=0;
	}

	void set_to_zero(ctype** ta) {
		int ic1,ic2;
		ctype* l01;
		for(ic1=0; ic1<n; ++ic1) {
			l01=ta[ic1];
			for(ic2=0; ic2<n; ++ic2) {
				l01[ic2]=0;
			}
		}
	}

	void set_to_identity(ctype** ta) {
		int ic1,ic2;
		ctype* l01;
		for(ic1=0; ic1<n; ++ic1) {
			l01=ta[ic1];
			for(ic2=0; ic2<n; ++ic2) {
				if(ic2==ic1) {
					l01[ic2]=1.;
				} else {
					l01[ic2]=0;
				}
			}
		}
	}

	void matrix_copy(ctype** lin1,ctype** lout) {
		int ic1,ic2;
		ctype* lin10;
		ctype* lout0;
		for(ic1=0; ic1<n; ++ic1) {
			lout0=lout[ic1];
			lin10=lin1[ic1];
			for(ic2=0; ic2<n; ++ic2) {
				lout0[ic2]=lin10[ic2];
			}
		}
	}

	void matrix_copy_a(ctype** lin1,ctype** lout) {
		int ic1,ic2;
		ctype* lin10;
		for(ic1=0; ic1<n; ++ic1) {
			lin10=lin1[ic1];
			for(ic2=0; ic2<n; ++ic2) {
				lout[ic2][ic1]=std::conj(lin10[ic2]);
			}
		}
	}

	void matrix_mult_nn(ctype** lin1,ctype** lin2,ctype** lout) {
		int ic1,ic2,ic3;
		ctype* lin10;
		ctype* lout0;
		ctype tout;
		for(ic1=0; ic1<n; ++ic1) {
			lout0=lout[ic1];
			lin10=lin1[ic1];
			for(ic2=0; ic2<n; ++ic2) {
				tout=0;
				for(ic3=0; ic3<n; ++ic3) {
					tout+=lin10[ic3]*lin2[ic3][ic2];
				}
				lout0[ic2]=tout;
			}
		}
	}

	void matrix_mult_nn(ctype** lin1,ctype** lin2,ftype rs,ctype** lout) {
		int ic1,ic2,ic3;
		ctype* lin10;
		ctype* lout0;
		ctype tout;
		for(ic1=0; ic1<n; ++ic1) {
			lout0=lout[ic1];
			lin10=lin1[ic1];
			for(ic2=0; ic2<n; ++ic2) {
				tout=0;
				for(ic3=0; ic3<n; ++ic3) {
					tout+=lin10[ic3]*lin2[ic3][ic2];
				}
				lout0[ic2]=rs*tout;
			}
		}
	}

	void matrix_mult_na(ctype** lin1,ctype** lin2,ctype** lout) {
		int ic1,ic2,ic3;
		ctype* lin10;
		ctype* lin20;
		ctype* lout0;
		ctype tout;
		for(ic1=0; ic1<n; ++ic1) {
			lout0=lout[ic1];
			lin10=lin1[ic1];
			for(ic2=0; ic2<n; ++ic2) {
				lin20=lin2[ic2];
				tout=0;
				for(ic3=0; ic3<n; ++ic3) {
					tout+=lin10[ic3]*std::conj(lin20[ic3]);
				}
				lout0[ic2]=tout;
			}
		}
	}

	void matrix_mult_na(ctype** lin1,ctype** lin2,ftype rs,ctype** lout) {
		int ic1,ic2,ic3;
		ctype* lin10;
		ctype* lin20;
		ctype* lout0;
		ctype tout;
		for(ic1=0; ic1<n; ++ic1) {
			lout0=lout[ic1];
			lin10=lin1[ic1];
			for(ic2=0; ic2<n; ++ic2) {
				lin20=lin2[ic2];
				tout=0;
				for(ic3=0; ic3<n; ++ic3) {
					tout+=lin10[ic3]*std::conj(lin20[ic3]);
				}
				lout0[ic2]=rs*tout;
			}
		}
	}

	void matrix_mult_an(ctype** lin1,ctype** lin2,ctype** lout) {
		int ic1,ic2,ic3;
		ctype* lout0;
		ctype tout;
		for(ic1=0; ic1<n; ++ic1) {
			lout0=lout[ic1];
			for(ic2=0; ic2<n; ++ic2) {
				tout=0;
				for(ic3=0; ic3<n; ++ic3) {
					tout+=std::conj(lin1[ic3][ic1])*lin2[ic3][ic2];
				}
				lout0[ic2]=tout;
			}
		}
	}

	void matrix_mult_an(ctype** lin1,ctype** lin2,ftype rs,ctype** lout) {
		int ic1,ic2,ic3;
		ctype* lout0;
		ctype tout;
		for(ic1=0; ic1<n; ++ic1) {
			lout0=lout[ic1];
			for(ic2=0; ic2<n; ++ic2) {
				tout=0;
				for(ic3=0; ic3<n; ++ic3) {
					tout+=std::conj(lin1[ic3][ic1])*lin2[ic3][ic2];
				}
				lout0[ic2]=rs*tout;
			}
		}
	}


	ctype get_trace(ctype** ta) {
		ctype res=0;
		for(int ic1=0; ic1<n; ++ic1) {
			res+=ta[ic1][ic1];
		}
		return res;
	}

	ctype get_det(ctype** ta) {
		ctype tdet=0;
		ctype ttdet=0;
		for(int i=0; i<nlcs; ++i) {
			ttdet=ta[0][lcs[i][0]]*(ftype)lcs[i][n];
			for(int j=1; j<n; ++j) {
				ttdet*=ta[j][lcs[i][j]];
			}
			tdet+=ttdet;
		}
		return tdet;
	}

	int n;
	esolve_h es;
	bessel bess;
	int** lcs;
	int nlcs;
	ftype eqprec;
	ctype** a;
	ftype* d;
	ftype* m;
};