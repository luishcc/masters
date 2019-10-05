#include <stdio.h>
#include <time.h>
#include <stdlib.h> 


int main(){
	clock_t t1,t2;
	double dt, sumdt;
	int i, j, num, cont, vez;
	float **a, *b, *x;
	FILE *fp;
	
	fp = fopen("coluna-c.txt", "w");
	num = 100;
	for(vez=1;vez<100;vez++){ /* mudar de 10 para 100 se for rodar ate 10^4 */
	
		b = (float*)malloc(num * sizeof(float));
		x = (float*)malloc(num * sizeof(float));
		a = (float**)malloc(num * sizeof(float*));
		for(i = 0; i < num; i++){
 			a[i] = (float*)malloc(num * sizeof(float));	
		}
 		

		for(cont=1;cont<10;cont++){
			sumdt = 0;
			
			for(i=1;i<num;i++){
				b[i]= 0.0;    
				x[i]= rand();
				for(j=1;j<num;j++){
					a[i][j] = rand();
				}
			}	
	
			t1 = clock();
			for(j=1;j<num;j++){
				for(i=1;i<num;i++){
					b[i] = b[i] + (a[i][j]*x[j]);
				}
			}
			t2 = clock();

			dt = ((double) (t2-t1))/ CLOCKS_PER_SEC;
			sumdt = sumdt + dt;

			printf("%d \n", cont);		
		}	
	
		sumdt = sumdt*0.1;

		printf("%8.8e ;  %i \n", sumdt, num);
		fprintf(fp, "%8.8e ;  %i \n", sumdt, num);

		num = num + 100;
	}

	fclose(fp);
}
