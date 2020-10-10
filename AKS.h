#ifndef __AKS_POL_H
#define __AKS_POL_H

#include<stdio.h>
#include<gmp.h>


typedef struct 
{
	mpz_t table;
	unsigned int size;
} __sieve;

typedef struct 
{
	mpz_ptr coef;
	unsigned int deg;
	unsigned char nil;
}	__mpz_poly;

typedef __mpz_poly mpX_t[1];
typedef __sieve sieve[1];


void 
mpX_init(mpX_t );
void 
mpX_clear(mpX_t );
void
mpX_set(mpX_t ,const mpX_t );
void 
mpX_set_ui(mpX_t  ,unsigned int );
void 
get_coef(mpz_t Z,const mpX_t , unsigned int );
unsigned long int 
get_coef_ui(const mpX_t , unsigned int );
void 
set_coef(mpX_t , mpz_t , unsigned int );
void 
set_coef_ui(mpX_t , unsigned int , unsigned int );
unsigned int 
mpX_isEqual(const mpX_t ,const mpX_t );
void 
mpX_print(const mpX_t );
void
mpX_mul(mpX_t ,const mpX_t ,const mpX_t );
void
mpX_add(mpX_t ,const mpX_t ,const mpX_t );
void 
mpX_mod_mul( mpX_t ,const mpX_t ,const mpX_t ,const mpz_t ,unsigned int );
void 
mpX_mod_pow(mpX_t , const mpX_t ,const mpz_t , const mpz_t,unsigned int );
void
initSieve(sieve );
void
clearSieve(sieve );
void
calcSieve(sieve ,unsigned int );
int
checkPrime(sieve ,unsigned int );
int 
AKS(mpz_t );
/* Initializing the polynomial */

void 
mpX_init(mpX_t pol)
{
	pol->coef = (mpz_ptr)NULL;
	pol->deg = 0;
	pol->nil = 0;
}

/* Clearing polynomial handles */

void 
mpX_clear(mpX_t pol)
{
	if(pol->nil==0)
		return;
	unsigned int d = pol->deg;
	int i;
	for(i=0;i<=d;i++)
		mpz_clear((pol->coef)+i);
	pol->coef = (mpz_ptr)NULL;
	pol->deg = 0;
	pol->nil = 0;
}

/* Copying one polynomial to another */

void
mpX_set(mpX_t A, const mpX_t B)
{
	mpX_clear(A);
	mpz_t c;
	mpz_init(c);

	unsigned int i,degB;
	degB = B->deg;
	for(i=0;i <= degB;i++)
	{
		get_coef(c,B,i);
		set_coef(A,c,i);
	}
	A->deg = degB;
	A->nil=1;
	mpz_clear(c);

}

/* Setting polynomial with one coefficient at 0'th place */

void 
mpX_set_ui(mpX_t pol ,unsigned int n)
{
	if(pol->nil == 1)
		mpX_clear(pol);
	
	pol->deg = 0;
	pol->nil = 1;
	pol->coef = (mpz_ptr)malloc(sizeof(mpz_t));
	if(pol->coef==NULL)
	{
		printf("Error in allocation\n");
		exit(1);
	}
	mpz_init_set_ui(pol->coef,n);
}

/* Getting i'th coefficient returns 0 if no coefficient exits */

void 
get_coef(mpz_t Z,const mpX_t pol, unsigned int i)
{
	unsigned int d;
	if(pol->nil == 0)
	{
		printf("Error in get_coef(...)\n");
		exit(1);
	}
	
	d = pol->deg;
	if(d < i)
	{
		mpz_set_ui(Z,0);
		return;
	}
	
	mpz_set(Z, (pol->coef)+i);
}

unsigned long int 
get_coef_ui(const mpX_t pol, unsigned int i)
{
	unsigned int d,coefficient;
	mpz_t Z;
	
	if(pol->nil == 0)
	{
		printf("Error in get_coef(...)\n");
		exit(1);
	}
	
	d = pol->deg;
	if(d < i)
		return 0UL;
	
	mpz_init_set(Z, (pol->coef)+i);
	coefficient = mpz_get_ui(Z);
	mpz_clear(Z);
	return coefficient;
}

/* Setting a coefficient at the i'th place */

void 
set_coef(mpX_t pol, mpz_t C, unsigned int k)
{
	int i;
	unsigned int d,size;
	if(mpz_cmp_ui(C,0)==0)
	{
		if(pol->nil==0)
			mpX_set_ui(pol,0);
		return;
	}
	
	d = pol->deg;
	if( pol->nil==1 && d >= k )
	{
		mpz_set((pol->coef)+k, C);
		pol->nil = 1;
		return;
	}
	else
	{
		size = k+1;
		pol->coef = (mpz_ptr)realloc(pol->coef,sizeof(mpz_t)*size);
		if(pol->coef==NULL)
		{
			printf("Error in allocation\n");
			exit(1);
		}
		for(i = d+1 ;i < k ; i++)
			mpz_init_set_ui((pol->coef)+i,0);
		
		mpz_init_set((pol->coef)+k, C);

		pol->deg=k;
		pol->nil = 1;
	}
}

void 
set_coef_ui(mpX_t pol, unsigned int C, unsigned int k)
{
	int i;
	unsigned int d,size;
	if(C==0)
	{
		if(pol->nil==0)
			mpX_set_ui(pol,0);
		return;
	}
	d = pol->deg;
	if( pol->nil==1 && d >= k )
	{
		
		mpz_init_set_ui((pol->coef)+k, C);
		
		
		
		pol->nil = 1;
		return;
	}
	else
	{
		size = k+1;
		pol->coef = (mpz_ptr)realloc(pol->coef,sizeof(mpz_t)*size);
		if(pol->coef==NULL)
		{
			printf("Error in allocation\n");
			exit(1);
		}
		for(i = d+1 ;i < k ; i++)
			mpz_init_set_ui((pol->coef)+i,0);
		mpz_init_set_ui((pol->coef)+k, C);
		pol->deg = k;
		pol->nil = 1;
	}
}

/* Checking if equal */

unsigned int 
mpX_isEqual(const mpX_t u,const mpX_t v)
{
	unsigned int i,degu,degv;
	degu = u->deg;
	degv = v->deg;

	if(degu != degv)
		return 0;
		
	mpz_t x,y;
	mpz_inits(x,y,NULL);

	for(i=0;i<=degu;i++)
	{
		get_coef(x,u,i);
		get_coef(y,v,i);
		if(mpz_cmp(x,y) != 0 )
		{
			mpz_clear(x);
			mpz_clear(y);
			return 0;
		}
	}
	mpz_clear(x);
	mpz_clear(y);
	return 1;
}

/* Printing a polynomial */

void 
mpX_print(const mpX_t pol)
{
	if(pol->nil == 0)
	{
		printf("\nEmpty polynomial\n\n");
		return;
	}
	unsigned int d = pol->deg;
	printf("degree : %d\n", d);
	for(int i=0;i<=d ;i++)
	{
		if(i==0)
		{
			gmp_printf("%Zd\n",pol->coef);
			continue;
		}
		if(mpz_cmp_ui(pol->coef + i,0) != 0)
			gmp_printf("%Zd X^%u\n",pol->coef + i, i);
	}
}

/* Multplication of two polynomial */

void
mpX_mul(mpX_t X, const mpX_t _u, const mpX_t _v)
{
	mpX_t u,v;
	
	mpX_init(u);
	mpX_init(v);
	mpX_set(u,_u);
	mpX_set(v,_v);
	mpX_clear(X);
	
	mpz_t x,y,sum;
	mpz_inits(x,y,sum,NULL);

	unsigned int i,j,degx, degu, degv;
	degu = u->deg;
	degv = v->deg;
	degx = degu + degv;

	for(i=0;i<=degx;i++)
	{
		mpz_set_ui(sum,0);

		for(j=0;j<=degu && j<=i;j++)
		{
			get_coef(x,u,j);
			get_coef(y,v,i-j);

			mpz_mul(x,x,y);
			mpz_add(sum,sum,x);
		}
		set_coef(X,sum,i);
	}
	X->deg = degx;
	X->nil=1;
	
	mpX_clear(u);
	mpX_clear(v);

	mpz_clear(x);
	mpz_clear(y);
	mpz_clear(sum);
}

/* Addition of two polynomial */

void
mpX_add(mpX_t X,const mpX_t _u,const mpX_t _v)
{
	mpX_t u,v;
	
	mpX_init(u);
	mpX_init(v);
	mpX_set(u,_u);
	mpX_set(v,_v);
	mpX_clear(X);
	
	mpz_t x,y;
	mpz_inits(x,y,NULL);

	unsigned int i,j,min,degx, degu, degv;
	degu = u->deg;
	degv = v->deg;
	degx = degu > degv ? degu : degv;
	min = degu < degv ? degu : degv;

	for(i=0;i<=min;i++)
	{
		get_coef(x,u,i);
		get_coef(y,v,i);
		mpz_add(x,x,y);
		set_coef(X,x,i);
	}
	while(i<=degu)
	{
		get_coef(x,u,i);
		set_coef(X,x,i);
		i++;
	}
	while(i<=degv)
	{
		get_coef(x,v,i);
		set_coef(X,x,i);
		i++;
	}
	X->deg = degx;
	X->nil=1;
	mpX_clear(u);
	mpX_clear(v);

	mpz_clear(x);
	mpz_clear(y);
}

/* Multiplication of two polynomial under (mod [X^r - 1]) and (mod n) */


/* [IMPORTANT : The assumption here is that the degree of _u and _v are at most (r-1)
                 since, that is all we require for this implementation.					]  */
void 
mpX_mod_mul( mpX_t X, const mpX_t _u, const mpX_t _v,const mpz_t n, unsigned int r)
{

	
	mpX_t u,v;
	mpz_t x,y,z,sum;
	unsigned int  i, j, k, maxdeg, degu, degv;
	
	mpX_init(u);
	mpX_init(v);

	mpX_set(u,_u);
	mpX_set(v,_v);
	
	mpX_clear(X);

	mpz_inits(x,y,z,sum,NULL);

	degu = u->deg;
	degv = v->deg;
	maxdeg = degu+degv;
	for(i=0 ; i < r ;i++)
	{
		mpz_set_ui(sum,0);

		for(j=0; j<=i ;j++ )
		{
			get_coef(x,u,j);
			if(mpz_cmp_ui(x,0) != 0)
			{
				get_coef(y,v,i-j);
				get_coef(z,v,i-j+r);
				
				mpz_add(y,y,z);
				mpz_mul(x,x,y);
				
				mpz_add(sum,sum,x);
			}
		}
		k = r+i;
		for(j=i+1 ; j<=k ; j++)
		{
			get_coef(x,u,j);
			if(mpz_cmp_ui(x,0) != 0)
			{
				get_coef(y,v,k-j);
				
				mpz_mul(x,x,y);
				
				mpz_add(sum,sum,x);
			}
		}
		mpz_mod(z,sum,n);
		set_coef(X,z,i);
		if(i >maxdeg && (mpz_cmp_ui(sum,0)==0))
			break;
	}
	X->nil=1;
	
	mpX_clear(u);
	mpX_clear(v);

	mpz_clear(x);
	mpz_clear(y);
	mpz_clear(z);
	mpz_clear(sum);

}

/* Making a polynomial to the power of some number */

void 
mpX_mod_pow(mpX_t R, const mpX_t _u,const mpz_t pow, const mpz_t n, unsigned int r)
{
	int bit,i;
	char * pattern;
	
	mpX_t u;
	mpX_init(u);
	mpX_set(u,_u);
	
	mpX_clear(R);
	mpX_init(R);
	mpX_set_ui(R,1);
	
	bit = mpz_sizeinbase(pow,2);
	i=0;
	pattern = mpz_get_str(NULL,2,pow);
	
	while( i<bit )
	{
		mpX_mod_mul(R,R,R,n,r);
		if(pattern[i]=='1')
			mpX_mod_mul(R,u,R,n,r);
			
		i++;
	}
	
	mpX_clear(u);	
}

/* Initialize sieve */

void
initSieve(sieve S)
{
	mpz_init_set_ui(S->table,0);
	S->size =2;
}

/* Clear Sieve */

void
clearSieve(sieve S)
{
	mpz_clear(S->table);
	S->size=0;
}

/* calculating sieve table */

void
calcSieve(sieve s,unsigned int r)
{
	unsigned int i,j,_size,size = s->size;
	if(size >= r)
		return;
	
	while(size < r)
	{
		_size = size;
		size *= 2;
		for(i=2; i<=_size ;i++)
		{
			if(!mpz_tstbit(s->table,i))
			{
				for(j=2*i ;j<=size; j+=i)
					mpz_setbit(s->table,j);
			}
		} 
	}
	s->size = size;
}

/* checking for small primes using seive prime */

int
checkPrime(sieve s,unsigned int r)
{
	if(s->size < r)
		calcSieve(s,r);
	return !mpz_tstbit(s->table,r);
}

/* AKS algorithm */

int 
AKS(mpz_t n)
{
	if(mpz_perfect_power_p(n))
		return 0;

	mpz_t x,_r,_sqrt,_lim;
	mpX_t lpoly,rpoly,__pol,__pol1,__pol2;

	mpz_inits(x,_r,_sqrt,_lim,NULL);
	
	sieve s;
	initSieve(s);

	unsigned int a,pow,r=2,logn,limit,res;
	logn = mpz_sizeinbase(n,2);
	limit = logn*logn;

	

	while(mpz_cmp_ui(n,r) > 0)
	{
		if(mpz_divisible_ui_p(n,r))
			return 0;
		int check = 0;


		if(checkPrime(s,r))
		{
			mpz_set_ui(_r,r);
			for(pow = 1; pow <= limit ;pow++)
			{

				mpz_powm_ui(x,n,pow,_r);

				if(mpz_cmp_ui(x,1)==0)
				{
					check = 1;
					break;
				}

			}
			if(check==0)
				break;

		}
		r++;
	}
	if(mpz_cmp_ui(n,r)==0)
		return 1;

	
	clearSieve(s);	
	
	mpz_set_ui(_r,r);
	mpz_sqrt(_sqrt,_r);
	mpz_add_ui(_sqrt,_sqrt,1);
	mpz_mul_ui(_lim,_sqrt,logn);

	mpX_init(lpoly);
	mpX_init(rpoly);
	mpX_init(__pol1);
	mpX_init(__pol2);
	mpX_init(__pol);

	mpz_mod(x,n,_r);
	set_coef_ui(__pol,1,mpz_get_ui(x));
	a=1;
	while(mpz_cmp_ui(_lim,a) > 0 )
	{
	
		#ifdef __VERBOSE_AKS
			if(a==1 || a%10==0)
			{
				printf(" %3.2lf%% \r",100*(float)a/mpz_get_ui(_lim));
				fflush(stdout);
			}
		#endif
		
		set_coef_ui(lpoly,1,1);
		set_coef_ui(lpoly,a,0);
		mpX_mod_pow(__pol1,lpoly,n,n,r);
		mpz_set_ui(x,a);
		mpz_mod(x,x,n);
		
		set_coef(rpoly,x,0);
		mpX_add(__pol2,rpoly,__pol);
		
		res = mpX_isEqual(__pol1,__pol2);

		if(res != 1)
		{
			mpX_clear(lpoly);
			mpX_clear(rpoly);
			mpX_clear(__pol);
			mpX_clear(__pol1);
			mpX_clear(__pol2);
			mpz_clear(x);
			mpz_clear(_r);
			mpz_clear(_sqrt);
			mpz_clear(_lim);
			return 0;
		}
		a++;
		
	}
	
	

	mpX_clear(lpoly);
	mpX_clear(rpoly);
	mpX_clear(__pol);
	mpX_clear(__pol1);
	mpX_clear(__pol2);
	mpz_clear(x);
	mpz_clear(_r);
	mpz_clear(_sqrt);
	mpz_clear(_lim);
	return 1;
}

#endif
