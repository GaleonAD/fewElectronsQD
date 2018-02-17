/* Driver for routine jacobi */

#include <stdio.h>
#include <math.h>
#define NRANSI
//#include "nr.h"
#include "nrutil.h"
#include "nrutil.c"
#include "jacobi.c"
#include "eigsrt.c"

void copyM( float ** M, float ** Mout, int N );
void tranM( float ** M, float ** Mout, int N );
void printM( float ** M, int N );
void printV( float * V, int N );
void multiplyM( float ** A, float ** B, float ** C,  int N );
float maxV( float * V, int N );
float maxMdifM( float ** M1, float ** M2, int Nx, int Ny );
float maxM( float ** M, int Nx, int Ny );


void copyV( float * M, float * Mout, int N );
float integral2D( float ** G1 , float ** G2, float dx, int Nx, float dy, int Ny );

///////////////////////////////////////////////////////////////////////////////////
//PARAMETRY I DEKLARACJE POMOCNICZE////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
// obsługa pliku, jednostki, parametry pudła obliczeniowego, studni
// 
///////////////////////////////////////////////////////////////////////////////////

int main(void)
{
	FILE * pFile;
	pFile = fopen("wynik","w");

// zmienne pomocnicze
	float x,y,x0,y0;	
	float Q,norma,r ;
	float HGj;
	int nrot; 
	int W = 3; 		// ilosc zapisywanych wektorów i wartości
	int l = 0; 		// na razie tylko stany l = 0
// jednostki
	double a_ = 18.8972; // długośc z nm na atomowe
	double e_ = 0.03649; // energia z eV na atomowe


// siatka i potencjal kropki
	double dx = 0.04*a_; 	// wielkosc kroku przestrzennego
	double Lx = 10.*a_;	// wilkosc pudla obliczeniowego (w nm)
	double dy = dx; 	// wielkosc kroku przestrzennego
	double Ly = Lx;		// wilkosc pudla obliczeniowego (w nm)
	double V0 = -25.*e_;	// potencjal w eV wewnatrz, zewnatrz jest 0
	double R =1.*a_;	// promien potencjalu (w nm)
	int Nx = 250;	//Lx/dx <- ale lepiej samemu wpisać :) punkty podzialu siatki
	int Ny = Nx; 	//Ly/dy;
	
////////////////////////////////////////////////////////////////////////////////////
//POTENCJAL KROPKI//////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
// ustalamy na skok
////////////////////////////////////////////////////////////////////////////////////

	float ** U = matrix(1,Nx,1,Ny);
	for (int i = 1; i <= Nx; ++i){ 
		for (int j = 1; j <= Ny; ++j) {
			if( fabs((i*dx-Lx/2.)*(i*dx-Lx/2.) 
				+ (j*dy-Ly/2.)*(j*dy-Ly/2.)) < R*R){
				U[i][j] = V0;
			}
			else{
				U[i][j]=0.;
			}	
		}
	}
////////////////////////////////////////////////////////////////////////////////////
//BAZA GAUSSIANOW///////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
// 
////////////////////////////////////////////////////////////////////////////////////
	printf("Tworzenie bazy...\n");	
	
	int Nbok = 7; // ilosc gaussianów na jeden bok
	int N = Nbok*Nbok;	// wielkosc tablicy gaussianów
	// tablica na baze
	float ** G[N];
	for( int i = 0; i < N ; ++i){
	 	G[i] = matrix(1,Nx,1,Ny);
	}
	// wypełnianie bazy <- MAGIA 
	//
	float sigma = 0.8; // 66% szer gaussianu <- trzeba zmieniać przy zmianie Nbok!!!
//	float bbb = 1.15*(Nbok-1.);
//	sigma = (-bbb+sqrt(bbb*bbb+32.*Lx))/16.;
	sigma = 0.078*Lx;
//	printf("sigma : %3f\n",sigma);
	x0 = 3*sigma; // odstep od krawedzi siatki <- trzeba zmieniac przy zmianie Nbok!!!
//	x0 = 1.5*Lx/10.; // odstep od krawedzi siatki <- trzeba zmieniac przy zmianie Nbok!!!
	y0 = x0; 		// j w
	float ddx = (Lx-2.*x0)/(Nbok-1); // odstep miedzy gaussianami
	float ddy = (Lx-2.*y0)/(Nbok-1); // jw
////// tu powinno być około 1.15!!!!!!!!!!!!!!!
	//printf("sigma : %3f\n",ddx);
	printf("sigma : %3f\n",ddx/sigma);
/////  wypełnianie bazy!!!
	for( int kx = 0; kx < Nbok; kx++){
	for( int ky = 0; ky < Nbok; ky++){
		for (int i = 1; i <= Nx; ++i){ 
			for (int j = 1; j <= Ny; ++j) {
				x = i*dx - x0 - kx*ddx;
				y = j*dy - y0 - ky*ddy;
				G[kx+Nbok*ky][i][j] = exp( (-x*x -y*y)/(sigma*sigma) );
			}
		}
	}}

////////////////////////////////////////////////////////////////////////////////////
//MACIERZ S/////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
// macierz przekrywania i jej diagonalizacje liczymy jeden raz, gdyz ona nie zalezy
// od potencjalu wejsciowego, czyli w peli sie nie zmienia
////////////////////////////////////////////////////////////////////////////////////
	printf("Tworzenie macierzy przekrywania...\n");	

        float ** S = matrix(1,N,1,N);
// wypełnianie macierzy przekrywania 
        for (int i = 1; i <= N; ++i){
                for (int j = 1; j <= N; ++j){
                       S[i][j] = integral2D( G[i-1], G[j-1], dx , Nx, dy, Ny );
                }
        }

// diagonalizacja S
        int n_rot_S; // do procedury jacobi
        float ** A = matrix(1,N,1,N);
        float ** AT = matrix(1,N,1,N);
        float  *eigenVS = vector(1,N); // wektor wartości własnych S
        
        // diagonalizacja S <- dostajemy macierz A przejścia 
        jacobi(S, N, eigenVS, A, &n_rot_S);
        // normalizacja A, tzn sprowadzenie macierzy A 
        // do postaci dającej iloczyn ATSA=I
        for (int i = 1; i <=N; ++i){
                for( int j = 1; j <= N; ++j){
                        A[i][j] = A[i][j] / sqrt(eigenVS[j]);
                }
        }
        // transpozycja A
        tranM( A, AT, N);
////////////////////////////////////////////////////////////////////////////////////
//PETLA OBLICZEN////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
// petla w ktorej relaksuje sie potencjal hartree
// deklaracja macierzy hamiltonianu i innych pomocniczych
////////////////////////////////////////////////////////////////////////////////////
	printf("Start pętli...\n");	
	// potencjal
	float alpha = 1e-2; // mały dodatek potencjału Hartree
	float ** Vin = matrix(1,Nx,1,Ny); // potencjał wejściowy (bez U kropki)
	float ** VH = matrix(1,Nx,1,Ny); // pomocnicze
	float ** Vout = matrix(1,Nx,1,Ny); // potencjał uzyskiwany z koncentracji
	for( int i = 1; i <= Nx; ++i){
	for( int j = 1; j <= Ny; ++j){
		Vin[i][j] = VH[i][j] = Vout[i][j] = 0.; // startowe zerowanie
	}}
	// do macierzy hamiltonianu
	float ** HA = matrix(1,N,1,N); // macierz tymczasowa pomocnicza
	int n_rot_H = 0; // do procedury jacobi
	float ** H = matrix(1,N,1,N); 	 // macierz hamiltonianu
	float ** ATHA = matrix(1,N,1,N); // na macierz hamiltonianu przygotowana do diagonalizacji
	float  *lambda = vector(1,N); // wektor wartości własnych H
	float ** B = matrix(1,N,1,N);	// macierz wektorów własnych ATHA
	float ** C = matrix(1,N,1,N);	// macierz wektorów własnych ATHA
	// do funkcji falowej 
	float ** chi0 = matrix(1,Nx,1,Ny); // zerowa funkcja falowa
	float ** chi1 = matrix(1,Nx,1,Ny); // 1. f. f.
	float ** chi2 = matrix(1,Nx,1,Ny); // 2. f. f.
	float epsNorma = 1e-3; // sprawdzana dokładność do normalizacji funkcji
	// do koncentracji i r. poissona
	float ** ne = matrix(1,Nx,1,Ny); // koncentracja elektornowa
	int Eiter; // do automatyzacji sumowania różnych funkcji falowych 
	int relaksMAX = 1e3; // liczba iteracji w metodzie relaksacji
	float epsVP = 1e-3; // sprawdza, czy metoda relaksacji działa wystarczająco dobrze
	// petla
	float epsVH = 1e-3; // = VHold-VHnew = 0. czyli warunek zaończenia pętli
	int iter = 1.; // iterator
	int iter_MAX = 1; // maksymalna ilosć iteracji 

	// start petli

	do{	
////////////////////////////////////////////////////////////////////////////////////
//POTENCJAL/////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
// dodawanie potencjalu hartree
////////////////////////////////////////////////////////////////////////////////////
	// liczenie nowego potencjalu wejściowego hartree
	//copyV(VHnew, VHold, N); 
	for( int i = 1; i <= Nx; ++i){
	for( int j = 1; j <= Ny; ++j){
		Vin[i][j] =  (1.-alpha)*Vin[i][j] + alpha*Vout[i][j];
	}}

////////////////////////////////////////////////////////////////////////////////////
//MACIERZ HAMILTONINANU/////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
// 
////////////////////////////////////////////////////////////////////////////////////
	// wypełmianie macierzy hamiltonianu H
		for (int i = 1; i <= N; ++i){
			for (int j = 1; j <= N; ++j){
				H[i][j] = 0.; // zerowanie poprzedniej iteracji
				// calka
				for (int xi = 2; xi < Nx; ++xi){
				for (int yj = 2; yj < Ny; ++yj){
					HGj = 0.;
					HGj += -( G[j-1][xi-1][yj] + G[j-1][xi+1][yj] + G[j-1][xi][yj-1] + G[j-1][xi][yj+1]
						- 4.*G[j-1][xi][yj])/(2.*dx*dx); // operator pedu
					HGj += U[xi][yj]*G[j-1][xi][yj]; //+ potencjal U kropki
					HGj += Vin[xi][yj]*G[j-1][xi][yj]; //+ potencjal Hartree
					H[i][j] += G[i-1][xi][yj]*HGj*dx*dy; // całkowanie 
				}} 
			}
		}
	
	// tworzenie macierzy ATHA 
		multiplyM( H, A, HA, N); 	// mnozenie H*A = HA
		multiplyM( AT, HA, ATHA, N);  // mnozenie AT*HA
		
	// diagonalizacja ATHA 
		jacobi(ATHA, N, lambda, B, &n_rot_H);
        	eigsrt(lambda,B,N);  //sortowanie energii od największej do najmniejszej
		
////////////////////////////////////////////////////////////////////////////////////
//FUNKCJE FALOWE ///////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
// obliczanie kilu funkcji własnych
////////////////////////////////////////////////////////////////////////////////////
	
		multiplyM( A, B, C, N); // macierz C= AB <- zawiera współczynniki 
					// do funkcji falowych
		//funkcje falowe
		for( int i = 1; i <= Nx; ++i){
		for( int j = 1; j <= Ny; ++j){
			chi0[i][j] = 0.; // zerowanie poprzedniej iteracji
			for( int k = 1; k<=N; ++k ){
				chi0[i][j] +=  G[k-1][i][j]* C[k][N];
			}
		}}
		for( int i = 1; i <= Nx; ++i){
		for( int j = 1; j <= Ny; ++j){
			chi1[i][j] = 0.; // zerowanie poprzedniej iteracji
			for( int k = 1; k<=N; ++k ){
				chi1[i][j] +=  G[k-1][i][j]* C[k][N-1];
			}
		}}
		for( int i = 1; i <= Nx; ++i){
		for( int j = 1; j <= Ny; ++j){
			chi2[i][j] = 0.; // zerowanie poprzedniej iteracji
			for( int k = 1; k<=N; ++k ){
				chi2[i][j] +=  G[k-1][i][j]* C[k][N-2];
			}
		}}

		// normalizacja dla pewności
		norma = integral2D( chi0, chi0, dx, Nx, dy, Ny);
		norma = sqrt(norma);
		// jeśli norma nie jest zachowana to normalizuje 
		if( fabs(norma-1.) < epsNorma ){
			for( int i = 1; i <= Nx; ++i){
			for( int j = 1; j <= Ny; ++j){
				chi0[i][j] /= norma;
			}}
		}
		norma = integral2D( chi1, chi1, dx, Nx, dy, Ny);
		norma = sqrt(norma);
		// jeśli norma nie jest zachowana to normalizuje 
		if( fabs(norma-1.) < epsNorma ){
			for( int i = 1; i <= Nx; ++i){
			for( int j = 1; j <= Ny; ++j){
				chi1[i][j] /= norma;
			}}
		}
		norma = integral2D( chi2, chi2, dx, Nx, dy, Ny);
		norma = sqrt(norma);
		// jeśli norma nie jest zachowana to normalizuje 
		if( fabs(norma-1.) < epsNorma ){
			for( int i = 1; i <= Nx; ++i){
			for( int j = 1; j <= Ny; ++j){
				chi2[i][j] /= norma;
			}}
		}

////////////////////////////////////////////////////////////////////////////////////
//KONCENTRACJA//////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
// obliczanie koncentracji elektronowej
////////////////////////////////////////////////////////////////////////////////////
	
		// zerowanie poprzedniej koncentracji
		Eiter = 0; // licznik uwzględnionych funkcji falowych = ładunek!!
		for( int i = 1; i <= Nx; ++i){
		for( int j = 1; j <= Ny; ++j){
			ne[i][j] = 0.;	
		}}
		// stan (11)  => dwie mozliwosci spinu
		for( int i = 1; i <= Nx; ++i){
		for( int j = 1; j <= Ny; ++j){
			ne[i][j] += 2.*chi0[i][j]*chi0[i][j];	
		}}
		Eiter += 2.;
		// stan (21)  => dwie mozliwosci spinu
		for( int i = 1; i <= Nx; ++i){
		for( int j = 1; j <= Ny; ++j){
			ne[i][j] += 2.*chi1[i][j]*chi1[i][j];	
		}}
		Eiter += 2.;
		// stan (12)  => dwie mozliwosci spinu
	//	for( int i = 1; i <= Nx; ++i){
	//	for( int j = 1; j <= Ny; ++j){
	//		ne[i][j] += 2.*chi1[i][j]*chi1[i][j];	
	//	}}
	//	Eiter += 2.;

////////////////////////////////////////////////////////////////////////////////////
//POISSON///////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
// obliczanie potencjału Hertree
////////////////////////////////////////////////////////////////////////////////////
// pamietamy e=eps*4*pi = 1
// normalna metoda z krokiem Eulera

		//  metoda relaksacji GS, zerowe warunki brzegowe
		for( int i = 1; i <= Nx; ++i){
		for( int j = 1; j <= Ny; ++j){
			VH[i][j] = 0.;
		}}

		for( int k = 1; k<=relaksMAX; ++k){
			for( int i = 2; i < Nx; ++i){
			for( int j = 2; j < Ny; ++j){
				VH[i][j] = (VH[i-1][j]+VH[i+1][j]
					+VH[i][j-1]+VH[i][j+1]+ne[i][j]*dx*dy)/4.;
			}}
		}

	
	// jedna dodatkowa pętla w celu upewnienia się co do dokłądnośći metody relaksacji
		for( int i = 1; i <= Nx; ++i){
		for( int j = 1; j <= Ny; ++j){
			Vout[i][j] = VH[i][j];
		}}

			for( int i = 2; i < Nx; ++i){
			for( int j = 2; j < Ny; ++j){
				VH[i][j] = (VH[i-1][j]+VH[i+1][j]
					+VH[i][j-1]+VH[i][j+1]+ne[i][j]*dx*dy)/4.;
			}}

	
		// komunikat gdy metoda relaksacji nie obliczyła odpowiedniej dokładności
		// brana dokładność względna
		// wyszuje różnicę między VH i Vout i dzieli przez max VH
		if(maxMdifM(VH, Vout, Nx, Ny)/maxM(VH, Nx, Ny) > epsVP ){
			printf("Uwaga! Metoda relaksacji nie zbiega!!");
		}
		


	
////////////////////////////////////////////////////////////////////////////////////
// Obliczanie potencjału wymiany, na razie nie//////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////////
		

////////////////////////////////////////////////////////////////////////////////////
//KONIEC - PETLI OBLICZEN///////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
// petla w ktorej relaksuje sie potencjal hartree
// sprawdzam maksymalna roznice
////////////////////////////////////////////////////////////////////////////////////
		// | różnica V |
		// stop petli po zbyt wielu iteracjach  
		if( iter >= iter_MAX ){
			printf("To many iterations of procedure, broken after %d ",iter);
			printf("iterations when VHold-VHnew = %f\n",maxMdifM(Vin, Vout, Nx, Ny));
			break;
		}// jeśli nie stop to zwiększenie iteratora
		else ++iter;
		
		// takie tam dla ułatwienia oszacowania czasu
		if( iter == iter_MAX/4 ){
			printf("%d/%d...\n",iter,iter_MAX);}
		if( iter == iter_MAX/2 ){
			printf("%d/%d...\n",iter,iter_MAX);}
		if( iter == 3*iter_MAX/4 ){
			printf("%d/%d...\n",iter,iter_MAX);}
		
		//koniec pętli jeśli VHartree się już zmienia 
	}while ( maxMdifM(Vin, Vout, Nx, Ny) > epsVH );

	// wypisywanie ilości wykonanych iteracji
	printf("iterations = %d", iter);
	printV(lambda, N); 	// wektor końcowych wartości własnych 
				// pasuje aby 3 ostatnie wartości własne 
				// były ujemne (związane w studni!!
				// można zastąpić wypisywaniem ostatnich trzech
	printf("\n\n");

	// zapisywanie wyników!!
	for (int i=1;i<=Nx;i++){ 
		for (int j=1;j<=Ny;j++){
		fprintf(pFile,"%f\t ", i*dx);
		fprintf(pFile,"%f\t ", j*dx);
		fprintf(pFile,"%f\t ", U[i][j]);
		fprintf(pFile,"%f\t ", Vin[i][j]);
		fprintf(pFile,"%f\t ", chi0[i][j]);
		fprintf(pFile,"%f\t ", chi1[i][j]);
		fprintf(pFile,"%f\t ", chi2[i][j]);
		for( int k = 0; k < N; ++k){
		//	fprintf(pFile,"%f\t ", G[k][i][j]);
		}
		fprintf(pFile,"\n");
		}
	}
	

// zamykanie 
	fclose(pFile);
	for( int k = 0; k < N; ++k){
		free_matrix(G[k],1,Nx,1,Ny);
	}	
	free_matrix(S,1,N,1,N);
	free_matrix(A,1,N,1,N);
	free_matrix(AT,1,N,1,N);
	free_vector(eigenVS,1,N);
	free_matrix(Vin,1,Nx,1,Ny);
	free_matrix(Vout,1,Nx,1,Ny);
	free_matrix(VH,1,Nx,1,Ny);
	free_matrix(H,1,N,1,N);
	free_matrix(HA,1,N,1,N);
	free_matrix(ATHA,1,N,1,N);
	free_matrix(C,1,N,1,N);
	free_matrix(B,1,N,1,N);
	free_vector(lambda,1,N);
	free_matrix(chi0,1,Nx,1,Ny);
	free_matrix(chi1,1,Nx,1,Ny);
	free_matrix(chi2,1,Nx,1,Ny);
	free_matrix(ne,1,Nx,1,Ny);
	return 0;
}


////////////////////////////////////////////////////////////////////////////////////
//FUNKCJE///////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
// kopiowanie, transpozycja, wypisywanie, mnorzenie macierzy, itd
////////////////////////////////////////////////////////////////////////////////////
// kopiowanie macierzy
void copyM( float ** M, float ** Mout, int N ){
        for ( int i = 1; i <=N; ++i){
                for ( int j = 1; j <=N; ++j){
                        Mout[i][j] = M[i][j];
                }
        }
}
// transponowanie macierz
void tranM( float ** M, float ** Mout, int N ){
        for ( int i = 1; i <=N; ++i){
                for ( int j = 1; j <=N; ++j){
                        Mout[j][i] = M[i][j];
                }
        }
}
// wypisanie macierzy
void printM( float ** M, int N ){
        for ( int i = 1; i <=N; ++i){
                printf("[ ");
                for ( int j = 1; j <=N; ++j){
                        printf("%f \t",M[i][j]);
                }
                printf(" ]\n");
        }
}
// kopiowanie wektora
void copyV( float * Vin, float * Vout, int N ){
        for ( int i = 1; i <=N; ++i){
                Vout[i] = Vin[i];
        }
}
// całkowanie iloczynu G1 i G2 po całej płaszczyźnie
float integral2D( float ** G1 , float ** G2, float dx, int Nx, float dy, int Ny ){
        float out = 0.;
	for ( int i = 1; i <=Nx; ++i){
        	for ( int j = 1; j <=Ny; ++j){
        		out += G1[i][j]*G2[i][j]*dx*dy;
		}
        }
	return out;
}
// wypisywanie wektora
void printV( float * V, int N ){
                printf("[ ");
                for ( int j = 1; j <=N; ++j){
                        printf("%f \t",V[j]);
                }
                printf(" ]\n");
}
// mnożenie macierzy C = AB
void multiplyM( float ** A, float ** B, float ** C,  int N ){
        float suma;
        for ( int i = 1; i <=N; ++i){
                for ( int j = 1; j <=N; ++j){
                        suma = 0;
                        for( int k = 1; k <= N; ++k){
                                suma += A[i][k]*B[k][j];
                        }
                        C[i][j] = suma;
                }
        }
}
// wyszukanie wartości maksymalnej w sesie wartości bezwzględnej 
float maxV( float * V,  int N ){
        float max = V[1];
        for ( int i = 2; i <=N; ++i){
                if( V[i] > max ) max = V[i];
        }
        return max;
}
// j. w. dla różnicy macierzy
float maxMdifM( float **M1, float ** M2,  int Nx, int Ny ){
        float max = 0.;
        for ( int i = 1; i <=Nx; ++i){
        for ( int j = 1; j <=Ny; ++j){
                if( fabs(M1[i][j]-M2[i][j]) > max ){
			max = fabs(M1[i][j]-M2[i][j]);
		}
        }}
        return max;
}
// j. w. dla jednej macierzy
float maxM( float **M, int Nx, int Ny ){
        float max = 0.;
        for ( int i = 1; i <=Nx; ++i){
        for ( int j = 1; j <=Ny; ++j){
                if( fabs(M[i][j]) > max ){
			max = fabs(M[i][j]);
		}
        }}
        return max;
}




#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software 7&X*. */
