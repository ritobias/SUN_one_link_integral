#pragma once
#include <algorithm>
#include <complex>
#include <limits>
#include <cmath>
#include <iostream>

typedef double ftype;
typedef std::complex<ftype> ctype;

static const ftype _fmin=std::numeric_limits<ftype>::min();
static const ftype _fmax=std::numeric_limits<ftype>::max();
static const ftype _fprec=std::numeric_limits<ftype>::epsilon();
static const ftype _lfprec=std::log(_fprec);

namespace levi_civita {
	// functions for creating lists of the non-zero entries of the Levi-Civita symbol
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
			lcs[i]=new int[n+1]();
		}

		int icnt=0;
		int* tarr;
		tarr=new int[n];
		for(int i=0; i<n; ++i) {
			tarr[i]=i;
		}
		permute(lcs,tarr,n,icnt);
		delete[] tarr;
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
	// class that acts as functor that computes the coefficients of the characteristic
	// polynomial of a complex input matrix ta[][] using the Faddeev-Leverrier algorithm.
	// returns optionally also the matrix invers of ta[][]
public:
	charpoly(): n(0),lb(0)	{
		// default constructor
		// matrix size n will have to be defined by calling set_n() before functor can be used
	}

	charpoly(int tn): n(tn) {
		// constructor with provided matrix size n
		lb=new ctype**[2];
		for(int i=0; i<2; ++i) {
			new_matrix(lb[i],0);
		}
	}

	void set_n(int tn) {
		// (re-)allocates the temporary memomry required to perform the Faddeev-Leverrier iteration
		// on tnxtn matrices.
		if(tn>0) {
			// provided matrix size is valid:
			if(tn!=n) {
				// do something only if the requested matrix size, tn, differs from currently used one, n :
				if(lb!=0) {
					// free memory that has previously been allocated for a different matrix size:
					for(int i=0; i<2; ++i) {
						delete_matrix(lb[i]);
					}
					delete[] lb;
				}
				// allocate the memory for the new matrix size tn:
				n=tn;
				lb=new ctype**[2];
				for(int i=0; i<2; ++i) {
					new_matrix(lb[i],0);
				}
			}
		} else {
			// invalid matrix size:
			if(lb!=0) {
				// if memory is allocated, free it:
				for(int i=0; i<2; ++i) {
					delete_matrix(lb[i]);
				}
				delete[] lb;
				lb=0;
			}
			n=0;
		}
	}

	~charpoly() {
		for(int i=0; i<2; ++i) {
			delete_matrix(lb[i]);
		}
		delete[] lb;
		n=0;
	}

	void operator()(ctype* lc,ctype** ta) {
		// input: complex matrix ta[][] 
		// output: polynomial coeffs. --> lc[]  
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
		
		for(i=1; i<n; ++i) {
			ilb=1-ilb;
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
		}
	}

	void operator()(ctype* lc,ctype** tia,ctype** ta) {
		// input: complex matrix ta[][] 
		// output: polynomial coeffs --> lc[] 
		//         matrix inverse --> tia[][]
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

		for(i=1; i<n; ++i) {
			ilb=1-ilb;
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
		}

		ttlc=lc[0];
		for(i1=0; i1<n; ++i1) {
			tlb1=tlb[i1];
			ta1=tia[i1];
			for(i2=0; i2<n; ++i2) {
				ta1[i2]=tlb1[i2]/ttlc;
			}
		}
	}

	void operator()(ftype* lc,ctype** ta) {
		// assuming real charactersitic polynomial coefficients (from hermitian ta[i][j])
		// input: hermitian matrix ta[][] 
		// output: real polynomial coeffs --> lc[] 
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

		for(i=1; i<n; ++i) {
			ilb=1-ilb;
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
		}
	}

	void operator()(ftype* lc,ctype** tia,ctype** ta) {
		// assuming real charactersitic polynomial coefficients (from hermitian ta[i][j])
		// input: hermitian matrix ta[][] 
		// output: real polynomial coeffs --> lc[] 
		//         matrix inverse --> tia[][]
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

		for(i=1; i<n; ++i) {
			ilb=1-ilb;
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
		}

		ttlc=lc[0];
		for(i1=0; i1<n; ++i1) {
			tlb1=tlb[i1];
			ta1=tia[i1];
			for(i2=0; i2<n; ++i2) {
				ta1[i2]=tlb1[i2]/ttlc;
			}
		}
	}


private:
	void new_matrix(ctype**& ta,int init=-1) {
		// allocates memory for a complex nxn matrix and sets ta[][] to point at it
		ta=new ctype*[n];
		if(init==-1) {
			// no initialization
			for(int ic1=0; ic1<n; ++ic1) {
				ta[ic1]=new ctype[n];
			}
		} else if(init==0) {
			// initialize to zero matrix
			for(int ic1=0; ic1<n; ++ic1) {
				ta[ic1]=new ctype[n]();
			}
		} else {
			// initialize to identity matrix
			for(int ic1=0; ic1<n; ++ic1) {
				ta[ic1]=new ctype[n]();
				ta[ic1][ic1]=1.0;
			}
		}
	}

	void delete_matrix(ctype**& ta) {
		//deallocate the memory pointed at by ta[][] and sets the pointer to zero
		for(int ic1=0; ic1<n; ++ic1) {
			delete[] ta[ic1];
		}
		delete[] ta;
		ta=0;
	}

	int n; //nxn matrix size
	ctype*** lb; //pointer to temporary storage lb[2][n][n] used in Faddeev-Leverrier iteration
};

class bessel {
	// class serving as functor for computing modified bessel functions of first kind
	// for integer orders n and real argument x. 
	// Algorithm adapted from Numerical Recipes 2nd ed., Sec. 6.6
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
		!     this subroutine calculates the modified bessel function of 1st kind
		!     for integer order n and any real x.
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
		!     this subroutine calculates the modified bessel function of 1st kind
		!     for integer orders [n,...,n+nn-1] for any real x. 
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
	// class implementing the inplicitly double shifted QR alorithm
	// for obtaining the eigenvalues of a hermitian matrix as functor
public:
	esolve_h(): n(0),d(0),a(0),s(0),w(0) {
		// default constructor
	}

	esolve_h(int tn): n(tn),d(0),a(0) {
		// constructor with provided matrix size n
		s=new ftype[n-1]();
		w=new ctype[n];
	}

	void set_n(int tn) {
		if(tn>0) {
			if(tn!=n) {
				d=0;
				a=0;
				if(s!=0) {
					delete[] s;
				}
				if(w!=0) {
					delete[] w;
				}
				n=tn;
				s=new ftype[n-1]();
				w=new ctype[n];
			}
		} else {
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
			n=0;
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

template<class T,class fT>
class lu_det {
	// class serving as functor to compute the determinant of a matrix T ta[][] 
	// using LU decomposition (adapted from https://github.com/CFT-HY/HILA resp. Numerical Recipes, 2nd ed. p. 47 ff))
public:
	lu_det(): n(0),a(0),vv(0) {
		// default constructor
		// will require a call to set_n(tn) to specify the size tnxtn of the matrices on which the
		// functor will operate, so that the required temporary memory can be allocated.
	}

	lu_det(int tn) {
		// construtor with provided matrix size tn
		if(tn>0) {
			// valid matrix size --> allocate temporary memory
			n=tn;
			a=new T*[n];
			vv=new fT[n];
			for(int i=0; i<n; ++i) {
				a[i]=new T[n]();
			}
		} else {
			// invalid matrix size --> fall back to default construtor
			n=0;
			a=0;
			vv=0;
		}
	}

	~lu_det() {
		if(a!=0) {
			for(int i=0; i<n; ++i) {
				delete[] a[i];
			}
			delete[] a;
			a=0;
		}
		if(vv!=0) {
			delete[] vv;
			vv=0;
		}
		n=0;
	}

	void set_n(int tn) {
		if(tn>0) {
			if(tn!=n) {
				if(a!=0) {
					for(int i=0; i<n; ++i) {
						delete[] a[i];
					}
					delete[] a;
					a=0;
				}
				n=tn;
				a=new T*[n];
				for(int i=0; i<n; ++i) {
					a[i]=new T[n]();
				}
				if(vv!=0) {
					delete[] vv;
					vv=0;
				}
				vv=new fT[n];
			}
		} else {
			if(a!=0) {
				for(int i=0; i<n; ++i) {
					delete[] a[i];
				}
				delete[] a;
				a=0;
			}
			if(vv!=0) {
				delete[] vv;
				vv=0;
			}
			n=0;
		}
	}

	T operator()(T** ta) {
		// returns determinant of the provided matrix ta[][] (ta remains unchanged)
		int imax=-1;
		T det=0;
		int i,j,k;
		T* a1;
		T* ta1;
		T csum,dum;
		fT big,tmp;
		fT d=1;
		for(i=0; i<n; ++i) {
			//loop over rows to get scaling information
			a1=a[i];
			ta1=ta[i];
			big=0;
			for(j=0; j<n; ++j) {
				a1[j]=ta1[j]; // copy elements of ta[][] to a[][]
				tmp=std::norm(a1[j]);
				if(tmp>big) {
					big=tmp;
				}
			}

			if(big==0) {
				//if full row is 0, then det must be 0
				return det;
			}

			//save scaling:
			vv[i]=1.0/std::sqrt(big);
		}

		
		for(j=0; j<n; ++j) {
			//loop over rows

			//build lower triangle:
			for(i=0; i<j; ++i) {
				a1=a[i];
				for(k=0; k<i; ++k) {
					a1[j]-=a1[k]*a[k][j];
				}
			}

			//search for the pivot and start upper triangle:
			big=0;
			for(i=j; i<n; ++i) {
				a1=a[i];
				csum=a1[j];
				for(k=0; k<j; ++k) {
					csum-=a1[k]*a[k][j];
				}
				a1[j]=csum;

				tmp=vv[i]*std::abs(csum);
				if(tmp>=big) {
					imax=i;
					big=tmp;
				}
			}
			ta1=a[j];
			//swap rows if neede:
			if(j!=imax) {
				a1=a[imax];
				ta1=a[j];
				for(k=0; k<n; ++k) {
					std::swap(a1[k],ta1[k]);
				}
				d=-d;
				vv[imax]=vv[j];
			}

			//if diag is zero now, then det must be zero:
			if(std::abs(ta1[j])==0) {
				return det;
			}

			//divide by the pivot:
			if(j!=n-1) {
				dum=1.0/ta1[j];
				for(i=j+1; i<n; ++i) {
					a[i][j]*=dum;
				}
			}
		}

		//compute det:
		det=d;
		for(j=0; j<n; ++j) {
			det*=a[j][j];
		}

		return det;
	}

	void operator()(T** ta,fT& tldet, T& tsdet) {
		//computes for matrix ta[][] :  tldet=log(abs(det(ta))); tsdet=phase(det(ta))
		// uses the same LU algorithm as T operator()(T** ta)
		int imax=-1;

		tsdet=0; // initialize det to be zero in case algorithm finds that ta[][] is singular
		tldet=0;
		int i,j,k;
		T* a1;
		T* ta1;
		T csum,dum;
		fT big,tmp;
		fT d=1;
		for(i=0; i<n; ++i) {
			//loop over rows to get scaling information
			a1=a[i];
			ta1=ta[i];
			big=0;
			for(j=0; j<n; ++j) {
				a1[j]=ta1[j]; // copy elements of ta[][] to a[][]
				tmp=std::norm(a1[j]);
				if(tmp>big) {
					big=tmp;
				}
			}

			if(big==0) {
				//if full row is 0, then det must be 0
				return;
			}

			//save scaling:
			vv[i]=1.0/std::sqrt(big);
		}

		//loop over rows:
		for(j=0; j<n; ++j) {

			//build lower triangle:
			for(i=0; i<j; ++i) {
				a1=a[i];
				for(k=0; k<i; ++k) {
					a1[j]-=a1[k]*a[k][j];
				}
			}

			//search for the pivot and start upper triangle:
			big=0;
			for(i=j; i<n; ++i) {
				a1=a[i];
				csum=a1[j];
				for(k=0; k<j; ++k) {
					csum-=a1[k]*a[k][j];
				}
				a1[j]=csum;

				tmp=vv[i]*std::abs(csum);
				if(tmp>=big) {
					imax=i;
					big=tmp;
				}
			}
			ta1=a[j];
			//swap rows if neede:
			if(j!=imax) {
				a1=a[imax];
				ta1=a[j];
				for(k=0; k<n; ++k) {
					std::swap(a1[k],ta1[k]);
				}
				d=-d;
				vv[imax]=vv[j];
			}

			//if diag is zero now, then det must be zero:
			if(std::abs(ta1[j])==0) {
				return;
			}

			//divide by the pivot:
			if(j!=n-1) {
				dum=1.0/ta1[j];
				for(i=j+1; i<n; ++i) {
					a[i][j]*=dum;
				}
			}
		}

		//compute tldet and tsdet
		tsdet=d;
		tldet=0;
		for(j=0; j<n; ++j) {
			dum=a[j][j];
			tmp=std::abs(dum);
			if(tmp>0) {
				tldet+=std::log(tmp);
				tsdet*=dum/tmp;
			} else {
				tsdet=0;
			}
		}
		return;
	}

private:
	int n; //size of nxn matrix to operate on
	T** a; //pointer to temporary nxn matrix on which LU decomposition will be performed
	       //in order to leave input matrix unchanged
	fT* vv; //temporary array of length n to hold row-scaling factors
};

class sun_integrator {
public:
	sun_integrator(): a(0),n(0),d(0),chp(),cdet(),fdet(),tbr(0),ltbr(0),nltbr(0),tbrsq(0),ltbrsq(0),cp(0),cps(0),cpl(0),cpc(0),rfsqrt(0),rf(0),trs(0),nhlmax(0),lmax(0),ltabl(0),ltmf(0),lnormf(0),mmax(0),tbrpowtab(0),pal(0),pals(0),pall(0),tal(0),talc(0),nl(0),lfl(0) {

	}

	sun_integrator(int tn,int td=4): n(tn),d(td),chp(tn),cdet(tn),fdet(tn),tbr(0),ltbr(0),nltbr(0),tbrsq(0),ltbrsq(0) {
		_init();
	}

	void init(int tn,int td=4) {
		n=tn;
		d=td;

		chp.set_n(n);
		cdet.set_n(n);
		fdet.set_n(n);

		_init();
	}

	void _init() {
		if(n>0&&d>1) {
			new_matrix(a,0);
			cp=new ftype[n+1]();
			cps=new int[n+1]();
			cpl=new ftype[n+1]();
			cpc=new ctype[n+1]();
			rfsqrt=(ftype)(2*(d-1));
			rf=rfsqrt*rfsqrt;

			trs=1.0/rf;

			nhlmax=3;

			lmax=100;
			ltabl=new ftype[lmax]();
			ftype tlf=0;
			ltabl[0]=tlf;
			for(int l=1; l<lmax; ++l) {
				tlf+=std::log((ftype)l);
				ltabl[l]=(ftype)n*tlf;
			}
			ltmf=0;
			for(int i=1; i<=n; ++i) {
				ltmf+=std::log((ftype)i);
			}

			lnormf=0;
			for(int i=1; i<n; ++i) {
				lnormf+=(n-i)*std::log((ftype)i);
			}

			mmax=100;
			lfl=new ftype**[mmax];
			for(int im=0; im<mmax; ++im) {
				lfl[im]=new ftype*[lmax];
				for(int il=0; il<lmax; ++il) {
					lfl[im][il]=new ftype[n]();
				}
			}

			tbrpowtab=new ftype*[n];
			for(int i=0; i<n; ++i) {
				tbrpowtab[i]=new ftype[n]();
			}

			pal=new ftype[n]();
			pals=new int[n]();
			pall=new ftype[n]();
			tal=new ftype**[lmax];
			for(int l=0; l<lmax; ++l) {
				tal[l]=new ftype*[n];
				for(int j=0; j<n; ++j) {
					tal[l][j]=new ftype[n]();
				}
			}
			talc=new ftype**[lmax];
			for(int l=0; l<lmax; ++l) {
				talc[l]=new ftype*[n];
				for(int j=0; j<n; ++j) {
					talc[l][j]=new ftype[n]();
				}
			}
			nl=new ftype[n]();
		}
	}

	void set_beta(ftype beta) {
		tbr=beta*rfsqrt/(2.0*(ftype)n); //combine beta/(2*n) with rescaling factor rfsqrt of input matrix
		nltbr=(ftype)n*std::log(tbr);
		tbrsq=tbr*tbr;
		ltbr=std::log(tbr);
		ltbrsq=2.0*ltbr;

		// set up the lookup table for the values lfl[m][l][i]=log(f(m,l,j)), where f(m,l,j)=tbr^(2*m-j)*j!/((m+l)!*(m-j)!)
		// note that f(m,l,j) contains an extra factor of tbr^(-j), this due to a balancing procedure to improve the conditioning
		// number of the matrices tal[l][][], who's determinants need to be computed. The balancing is achived by multiplying
		// the each element tal[l][j][i] by tbr^(i-j)
		ftype lf0=0;
		for(int im=0; im<n; ++im) {
			ftype lfl0=lf0;
			for(int il=0; il<lmax; ++il) {
				ftype lfj0=lfl0;
				for(int ij=im; ij>=0; --ij) {
					lfl[im][il][ij]=lfj0;
					lfj0+=ltbr-std::log((ftype)(im-(ij-1)));
				}
				lfl0+=std::log((ftype)(1+il))-std::log((ftype)(1+im+il));
			}
			lf0+=ltbr-std::log((ftype)(1+im));
		}


		lf0+=ltbr;

		for(int im=n; im<mmax; ++im) {
			ftype lfl0=lf0;
			for(int il=0; il<lmax; ++il) {
				ftype lfj0=lfl0;
				for(int ij=n-1; ij>=0; --ij) {
					lfl[im][il][ij]=lfj0;
					lfj0+=ltbr-std::log((ftype)(im-(ij-1)));
				}
				lfl0+=std::log((ftype)(1+il))-std::log((ftype)(1+im+il));
			}
			lf0+=ltbrsq-std::log((ftype)(1+im))-std::log((ftype)(2+im-n));
		}
		
		// normally each column k of tal[l][j][k] would come with a rescaling factor tbrsq^(-k); but due to the aforementioned
		// balancing procedure, the factor is only tbr^(-k).
		ftype powfc=1.0;
		for(int k=0; k<n; ++k) {
			for(int j=0; j<n; ++j) {
				tbrpowtab[j][k]=powfc;
			}
			powfc/=tbr;
		}
	}

	~sun_integrator() {
		delete_matrix(a);
		delete[] cp;
		delete[] cps;
		delete[] cpl;
		delete[] cpc;
		delete[] ltabl;
		for(int im=0; im<mmax; ++im) {
			for(int il=0; il<lmax; ++il) {
				delete[] lfl[im][il];
			}
			delete[] lfl[im];
		}
		delete[] lfl;

		for(int i=0; i<n; ++i) {
			delete[] tbrpowtab[i];
		}
		delete[] tbrpowtab;

		delete[] pal;
		delete[] pals;
		delete[] pall;

		for(int l=0; l<lmax; ++l) {
			for(int j=0; j<n; ++j) {
				delete[] tal[l][j];
			}
			delete[] tal[l];
		}
		delete[] tal;
		for(int l=0; l<lmax; ++l) {
			for(int j=0; j<n; ++j) {
				delete[] talc[l][j];
			}
			delete[] talc[l];
		}
		delete[] talc;
		delete[] nl;
	}

	ftype operator()(ctype** staple) {

		// compute the hermitian matrix a=stapled.staple^{\dagger}:
		matrix_mult_na(staple,staple,trs,a);

		// since we won't need the matrix powers of a, we compute the coefficients of the
		// characteristic polynomial of a with the 
		chp(cp,a);

		int tlmax=lmax;
		int im,ij,il;
		ftype lrpf=0;
		ftype tacp0=std::abs(cp[0]);
		if(tacp0>_fprec) {
			lrpf=nltbr+0.5*std::log(tacp0);
			for(il=1; il<lmax; ++il) {
				if((ftype)il*lrpf-ltabl[il]<_lfprec) {
					tlmax=il;
					break;
				}
			}
		} else {
			tlmax=1;
		}

		for(im=0; im<n; ++im) {
			pal[im]=0;
			pals[im]=0;
			pall[im]=0;
			for(ij=0; ij<=im; ++ij) {
				for(il=0; il<tlmax; ++il) {
					tal[il][ij][im]=std::exp(lfl[im][il][ij]);
					talc[il][ij][im]=0;
				}
			}
			for(ij=im+1; ij<n; ++ij) {
				for(il=0; il<tlmax; ++il) {
					tal[il][ij][im]=0;
					talc[il][ij][im]=0;
				}
			}
		}
		pal[n-1]=1.0;
		pals[n-1]=1;
		pall[n-1]=0;

		int tch;
		int k;
		int nhl=nhlmax;
		ftype** lfc;
		ftype afcl;
		ftype ch,cho;
		ftype chl,chol,ch2l;
		ftype tt,ttal,texp;
		int chs,chos,ch2s;

		for(im=0; im<=n; ++im) {
			ch=cp[im];
			if(cp[im]==0) {
				cps[im]=0;
				cpl[im]=0;
			} else {
				if(cp[im]>0) {
					cps[im]=1;
					cpl[im]=std::log(cp[im]);
				} else {
					cps[im]=-1;
					cpl[im]=std::log(-cp[im]);
				}
			}
		}

		for(im=n; im<mmax; ++im) {
			lfc=lfl[im];
			afcl=lfc[0][0];
			tch=0;
			//ch=-pal[n-1]*cp[0];
			chs=-pals[n-1]*cps[0];
			if(chs!=0) {
				chl=pall[n-1]+cpl[0];
			} else {
				chl=0;
			}

			//cho=pal[0];
			chos=pals[0];
			chol=pall[0];

			//pal[0]=ch;
			pals[0]=chs;
			pall[0]=chl;

			if(chs!=0) {
				for(ij=0; ij<n; ++ij) {
					for(il=0; il<tlmax; ++il) {
						ttal=tal[il][ij][0];
						texp=std::exp(chl+lfc[il][ij]);
						if(texp>_fprec*std::abs(ttal)) {
							//tal[il][ij][0]+=chs*std::exp(chl+lfc[il][ij]); (with Kahan summation)
							texp*=(ftype)chs;
							tt=ttal+texp;
							if(std::abs(ttal)>=std::abs(texp)) {
								talc[il][ij][0]+=(ttal-tt)+texp;
							} else {
								talc[il][ij][0]+=(texp-tt)+ttal;
							}
							tal[il][ij][0]=tt;
							tch=1;
						}
					}
				}
			}

			for(k=1; k<n; ++k) {
				//ch2=-pal[n-1]*cp[k];
				ch2s=-pals[n-1]*cps[k];
				ch2l=pall[n-1]+cpl[k];

				//ch=cho+ch2;
				if(chos*ch2s==0) {
					if(chos==0) {
						chs=ch2s;
						if(chs!=0) {
							chl=ch2l;
						} else {
							chl=0;
						}
					} else {
						chs=chos;
						if(chs!=0) {
							chl=chol;
						} else {
							chl=0;
						}
					}
				} else if(ch2s==-chos&&std::abs(chol-ch2l)<_fprec) {
					chs=0;
					chl=0;
				} else {
					if(chol>ch2l) {
						chs=chos;
						if(chs!=0) {
							chl=chol+std::log(1.0+(ftype)(ch2s*chos)*std::exp(ch2l-chol));
						} else {
							chl=0;
						}
					} else {
						chs=ch2s;
						if(chs!=0) {
							chl=ch2l+std::log(1.0+(ftype)(chos*ch2s)*std::exp(chol-ch2l));
						} else {
							chl=0;
						}
					}
				}

				//cho=pal[k];
				chos=pals[k];
				chol=pall[k];

				//pal[k]=ch;
				pals[k]=chs;
				pall[k]=chl;

				if(chs!=0) {
					for(ij=0; ij<n; ++ij) {
						for(il=0; il<tlmax; ++il) {
							//tal[il][ij][k]+=chs*std::exp(chl+lfc[il][ij]); (with Kahan summation)
							ttal=tal[il][ij][k];
							texp=std::exp(chl+lfc[il][ij]);
							if(texp>_fprec*std::abs(ttal)) {
								texp*=(ftype)chs;
								tt=ttal+texp;
								if(std::abs(ttal)>=std::abs(texp)) {
									talc[il][ij][k]+=(ttal-tt)+texp;
								} else {
									talc[il][ij][k]+=(texp-tt)+ttal;
								}
								tal[il][ij][k]=tt;
								tch=1;
							}
						}
					}
				}
			}
			if(tch==0) {
				--nhl;
			}

			if(nhl>=0) {
				if(tch==1&&nhl<nhlmax) {
					nhl=nhlmax;
				}
			} else {
				std::cout<<"niter: "<<im<<std::endl;
				break;
			}
		}

		ctype ttdet=cdet(staple);
		ftype tth=std::arg(ttdet);

		il=0;
		for(ij=0; ij<n; ++ij) {
			for(k=0; k<n; ++k) {
				tal[il][ij][k]+=talc[il][ij][k];
				tal[il][ij][k]*=tbrpowtab[k][ij];
			}
		}

		if(0) {
			// print tal[0][][] matrix:
			std::cout<<std::endl;
			for(ij=0; ij<n; ++ij) {
				for(k=0; k<n; ++k) {
					std::cout<<tal[il][ij][k]<<" ";
				}
				std::cout<<std::endl;
			}
			std::cout<<std::endl;
		}

		fdet(tal[il],chl,ch);
		cho=(ftype)il*lrpf-ltabl[il];
		ftype maxexp=chl+cho;
		ftype tsum=ch;
		ftype tmaxexp;
		for(il=1; il<tlmax; ++il) {
			for(ij=0; ij<n; ++ij) {
				for(k=0; k<n; ++k) {
					tal[il][ij][k]+=talc[il][ij][k];
					tal[il][ij][k]*=tbrpowtab[k][ij];
				}
			}
			
			fdet(tal[il],chl,ch);
			cho=(ftype)il*lrpf-ltabl[il];
			tmaxexp=chl+cho;
			if(tmaxexp>maxexp) {
				tsum*=std::exp(maxexp-tmaxexp);
				maxexp=tmaxexp;
			}
			tsum+=2.0*std::cos(il*tth)*ch*std::exp(tmaxexp-maxexp);
		}
		return std::log(tsum)+maxexp+lnormf;
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
				l01[ic2]=0;
			}
			l01[ic1]=1.0;
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

	ftype get_ldet_stable(ftype** ta) {
		int i,j,k;
		ftype* tv1;
		ftype* tv2;
		ftype ch,cho;
		ftype ldet=0;
		for(i=0; i<n; ++i) {
			tv1=ta[i];
			cho=get_norm(tv1);
			ldet+=std::log(cho);
			cho=1.0/cho;
			for(k=0; k<n; ++k) {
				tv1[k]*=cho;
			}
			for(j=i+1; j<n; ++j) {
				tv2=ta[j];
				ch=0;
				for(k=0; k<n; ++k) {
					ch+=tv1[k]*tv2[k];
				}
				for(k=0; k<n; ++k) {
					tv2[k]-=ch*tv1[k];
				}
			}
		}
		return ldet;
	}

	ftype get_norm(ftype* tv) {
		ftype res=0;
		ftype tres;
		for(int j=0; j<n; ++j) {
			tres=tv[j];
			res+=tres*tres;
		}
		return std::sqrt(res);
	}

	ftype get_normsq(ftype* tv) {
		ftype res=0;
		ftype tres;
		for(int j=0; j<n; ++j) {
			tres=tv[j];
			res+=tres*tres;
		}
		return res;
	}

	ctype get_trace(ctype** ta) {
		ctype res=0;
		for(int ic1=0; ic1<n; ++ic1) {
			res+=ta[ic1][ic1];
		}
		return res;
	}

	ctype get_det(ctype** ta) {
		chp(cpc,ta);
		if(n%2==0) {
			return cpc[0];
		} else {
			return -cpc[0];
		}
	}

	void get_inverse(ctype** ta,ctype** tia) {
		chp(cpc,tia,ta);
	}

	int n;
	int d;
	int lmax;
	int mmax;
	int nhlmax;
	ftype* ltabl;
	ftype** tbrpowtab;
	ftype*** lfl;
	ftype rf;
	ftype rfsqrt;
	ctype** a;
	ctype* cpc;
	ftype* cp;
	int* cps;
	ftype* cpl;
	ftype* pal;
	int* pals;
	ftype* pall;
	ftype* nl;
	ftype*** tal;
	ftype*** talc;
	ftype tbr;
	ftype ltbr;
	ftype nltbr;
	ftype trs;
	ftype tbrsq;
	ftype ltbrsq;
	ftype ltmf;
	ftype lnormf;
	charpoly chp;
	lu_det<ctype,ftype> cdet;
	lu_det<ftype,ftype> fdet;
};
