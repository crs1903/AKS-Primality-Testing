#include<stdio.h>
#include<stdlib.h>
#include<gmp.h>
#include<time.h>

#define __VERBOSE_AKS		/* define this before including AKS.h if want to see the progress */
#include"AKS.h"

int main()
{
	unsigned long int i,t;
	
	
	char str[1025],*bit;
	mpz_t k;
	mpz_init(k);
	
	printf("Enter a number : ");
	if(!scanf("%1024s",str))
	{
		printf("Error in Scaning");
		exit(1);
	}
	
	i=clock();
	mpz_set_str(k,str,10);
	
	gmp_printf("\nEntered number is : %Zd, bits : %d\n\n",k,mpz_sizeinbase(k,2));
	if(AKS(k)==1)
		printf("The number is prime\n");
	else
		printf("The number is not prime\n");
	t=clock();
	printf("\n\nTime taken ~ %lf sec\n",(double)(t-i)/CLOCKS_PER_SEC);
}
