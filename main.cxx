// *  Ludovico Mattavelli 1193955    Esame del 01/06/2022

// * Programma per il calcolo delle energie e degli autostati di una particella in una buca di potenziale infinita con un gradino di potenziale posto all'interno

// ** put this inside the file of the introduction
#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <stdlib.h>

using namespace std;

//definition of the variable
//const variables general
const int    I_case=1;                  // 0= funzione quadratica, 1=buca di potenziale con sin
const int    N_stati=3;                 // numero di autostati da calcolare
const int    Nh=15;                     // numero di parametri
const double V0=0.1;                    // altezza del gradino [Hartree]

//costanti per la definizione dell'asse x [Bohr]
const double x_min=0.;                  // valore minimo  x della buca di potenziale
const double x_max=10.;                 // valore massimo x della buca di potenziale
const double a=3.;                      // inizio del gradino
const double b=7.;                      // fine del gradino 
const int    Nx=200;                    // numero punti dell'asse x
const double hx=(x_max-x_min)/(Nx-1.);  // spaziatura tra i pti x

//Spaziatura nello spazio dei parametri
const double h_der=0.001;               // variazione per la derivata

// Soglie di tolleranza
const double toll1=0.0000001;              // tolleranza per ricerca di minimo lungo u * why? 
const double toll2=0.00000001;             // tolleranza per ricerca di minimo nel gradiente coniugato

// Incremento nel minimo lungo u
const double LL=0.01;                      


// Funzioni utilizzate
double  V(double);										//Potenziale a gradino
double  Norm(double f[Nx]);									//Normalizzazione degli autostati
double  Modulo_v(double h[Nh]);									//Modulo di un vettore
double  Energia(double B_h[Nh],int n_stato, double PSI[N_stati][Nx], double x[Nx]);			//Calcolo dell'energia
double  DerivParziale(double B_h[Nh], int r, int n_stato, double PSI[N_stati][Nx],double x[Nx]); 	//DerivataParziale



int main() 
{
	// definizione e inizializzazione delle variabili
	// variabili per i cicli for;
	int i=0;
	int j=0;
	int k=0;
	int n_stato=0;
	
	//Definisco i file di output
	ofstream fout[N_stati];                     
	string nome_file[N_stati]; 
	
	//apro il file che continene i nomi che assegnerò ai files contenenti gli autovettori
  	string inputfilename1;
	ifstream fin1;
	fin1.open("Nomi_files_autovettori.txt"); 

	if(!fin1)
	{
		cout << "Errore in apertura del file1 di input." << endl;   
	   	return -1;
	}

	string file[N_stati];
	for(int i=0;i<N_stati;i++)
	{
	   	fin1>>file[i];
	}
	
        //Array per i parametri:               
	double B[Nh]={0.};                        // Condizioni iniziali del metodo del gradiente coniugato  
   	double B1[Nh]={0.};	    			     // Interazione i+1 
   	double B0[Nh]={0.};				     // Inizio ciclo 1
   
	//Vettore asse x
	double x[Nx];                                  
	for(i=0;  i<Nx; i++)
    	{
      		x[i]=x_min+((double) i)*(x_max-x_min)/(Nx-1.);
     	}
     	
     	
     	//Variabili per il metodo del gradiente coniugato
     	//Gradienti e vettori 
	double g0[Nh],h0[Nh], h1[Nh], g1[Nh]; 	     // 0=interazione i_esima, 1=interazione i+1_esima
	double u[Nh];				     // Versore relativo ad h1
	//Altre variabili   	
   	int mmm=0;                                  // numero di iterazioni del metodo del gradiente coniugato
    	int nnn=0;                                  // numero di iterazioni nella ricerca lungo la linea u
    	double GammaNum=0.;
    	double GammaDen=0.;
    	double Gamma=0.;
    	double Modulo_h1=0.;
    	double NormB_Denom=0.;
    	double De=0.;				     // Derivata lungo h1
	
	
	// Variabili per il salvataggio dello stato fondamentale
	double PSI[N_stati][Nx]={0.};			     // Autofunzioni
	double B1_sf[N_stati][Nh]={0.};			     // Parametri delle autofunzioni
	double f[Nx]={0.};
	double Norm_f=0.;
	double Coeff[N_stati]={0.};
	
	//Array per salvare i dati
  	double Autovalori[N_stati];                    // Energia per ogni stato
	double Parametri_Finali[N_stati][Nh];       // Parametri per ogni stato di eccitazione
	
	double Ortog=0.;			     // Verifica dell'ortogonalità
	
	
// Ciclo 0: per i diversi stati di eccitazione
for(n_stato=0; n_stato<N_stati; n_stato++)
{
	// Inizializzo tutti i parametri a 0
	for(j=0; j<Nh; j++)
        {
        	B[j]=0.;
        }
	// Pongo il parametro n_stato pari a 1
	B[n_stato]=1.;

	cout<<"Valore iniziale Energia stato "<<n_stato<<" : "<<Energia(B,n_stato,PSI,x)<<endl;
	
	// reset numero di iterazioni
	mmm=0;      
	
    	// Ciclo 1: metodo del gradiente coniugato
	do
	{
		// Salvo i parametri prima del ciclo
		for(i=0; i<Nh; i++)
		{
			B0[i]=B[i];
		}
		// Calcolo g_i+1=-Grad_Energia
	 	for (j=0; j<Nh; j++)
        	{
        		g1[j]=-DerivParziale(B,j,n_stato,PSI,x);     
        	}
        	
	   	// Calcolo h_i+1	
		if(mmm<1)			// lo eseguo solo la prima volta
		{
			for(j=0; j<Nh; j++)
			{
				h1[j]=g1[j];
			}
		}
		
		else				// dalla seconda iterazione
		{		
			//Calcolo il fattore Gamma
			GammaNum=0.;
			GammaDen=0.;
			for(j=0; j<Nh; j++)
			{
				GammaNum+=(g1[j]-g0[j])*g1[j];
				GammaDen+=(g0[j]*g0[j]);
			}
			
            		
            		Gamma=GammaNum/GammaDen;
			for(j=0; j<Nh; j++)
			{
				h1[j]=g1[j]+h0[j]*Gamma;
			}
		}
		
		// Calcolo il versore u di h1
		Modulo_h1=Modulo_v(h1);
		for(j=0; j<Nh; j++)
		{
		   	u[j]=h1[j]/Modulo_h1;
		}
		
 		// reset numero di iterazioni
        	nnn=0;   
        	// ciclo 2: ricerca del minimo lungo la linea del versore u                   
 		do 
 		{
 			if(nnn<1)               	// la prima volta non faccio nulla
 			{
 				nnn++;
 			}
 				
 			else				// dalla seconda iterazione
 			{
 				//aggiorno i valori di B;
 				for(j=0; j<Nh; j++)
 				{
 					B[j]=B1[j];
 				}
 				
 				nnn++;  //Conto il numero di iterazioni del ciclo 2
 			}
 			
 			// calcolo il gradiente di Energia lungo u, come prodotto interno tra il gradiente e il versore u
            		De=0.;
 			for(j=0; j<Nh; j++)
 			{
 				De+=DerivParziale(B,j,n_stato,PSI,x)*u[j];
   			}
   			
 			for(j=0; j<Nh; j++)
 			{
   				B1[j]=B[j]-LL*De*u[j];  
   			}
   			
 		} while(abs(Energia(B,n_stato,PSI,x)-Energia(B1,n_stato,PSI,x))>toll1);   // fine del cilclo 2

		//Aggiorno i valori per il prossimo ciclo 1
		for(j=0; j<Nh; j++) 
		{
			g0[j]=g1[j];
			B[j] =B1[j];
			h0[j]=h1[j];
		}
		
		mmm++;		//Conto il numero di iterazioni del ciclo 2
		
		
	} while(abs(Energia(B0,n_stato,PSI,x)-Energia(B1,n_stato,PSI,x))>toll2);  // fine del ciclo 1
	//fine gradiente coniutgato
    
    	
    	if (I_case==1)  // normalizzo solo nel caso di sinusoidi
    	{
    		// normalizzo per il vettore B
    		NormB_Denom=1./Modulo_v(B1)*sqrt(2./x_max);
    		for(j=0; j<Nh; j++)
    		{
        		B1[j]=B1[j]*NormB_Denom;
    		}
    	}
    	
    	cout<<"Numero cicli effettuati "<<mmm<<endl;
        cout<<"Parametri finali "<<B1[0]<<" "<<B1[1]<<" "<<B1[2]<<" "<<B1[3]<<" "<<B1[4]<<" "<<endl;   
        cout<<"L'energia minima è: "<<Energia(B1,n_stato,PSI,x)<<endl;
	
    	// Salvo i parametri e il risultato finale
    	Autovalori[n_stato]=Energia(B1,n_stato,PSI,x);
    	for(j=0; j<Nh; j++)
    	{
    	    Parametri_Finali[n_stato][j]=B1[j];
    	}
    	
    	
   	//stampo l'autofunzione
	fout[n_stato].open(nome_file[n_stato],ios::out);
	fout[n_stato].precision(9);
	
   	f[Nx]={0.};
	Norm_f=0.;
	Coeff[N_stati]={0.};
	for (j=0; j<Nh; j++)
	{
	 	for(i=0;  i<Nx; i++)
		{
			f[i]+=B1[j]*sin((j+1)*x[i]*M_PI/x_max);
		}
    	}
    	//Normalizzo f
   	Norm_f=Norm(f);
    	for(i=0;i<Nx;i++)
    	{
    	    f[i]=f[i]/Norm_f;
    	}
    	
    	if(n_stato>0) 
    	{
    		// Calcolo del coeff
    		for(k=0; k<n_stato;k++)
    		{
			Coeff[k]=0.5*(f[0]*PSI[k][0]+f[Nx-1]*PSI[k][Nx-1])*hx;
			for(i=1;i<Nx-1;i++)
			{
		    		Coeff[k]+=f[i]*PSI[k][i]*hx;
	       	}
       	}
       	
       	for(k=0; k<n_stato;k++)
       	{
	       	for(i=0; i<Nx; i++)
	       	{
		   		f[i]=f[i]-PSI[k][i]*Coeff[k];
	       	}
	       	
	       	//Calcolo i parametri corretti
	       	for(j=0; j<Nh; j++)
			{
		    		B1[j]=B1[j]-Coeff[k]*B1_sf[k][j];
			}
        	}
        	
        	Norm_f=Norm(f);
        	for(i=0; i<Nx; i++)
       	{
           		f[i]=f[i]/Norm_f;
       	}
       	     	
    	}
    
    	//Salvo l'autostato per le successive interazioni
    	for(i=0; i<Nx; i++)
    	{
    	        PSI[n_stato][i]=f[i];
    	}
    	for(j=0; j<Nh; j++)
    	{
    	        	B1_sf[n_stato][j]=B1[j];
    	}
    
    	// Stampo sul file di output l'autostato
	for(i=0; i<Nx; i++)
	{
		fout[n_stato]<<x[i]<<" "<<f[i]<<endl;
	}
	fout[n_stato].close();


	//Verifica dell'ortogonalità
	if(n_stato==(N_stati-1))
	{	
		for(j=0;j<n_stato;j++)
		{
			for(k=0;k<n_stato;k++)
			{
				if(k!=j)
				{
					// Calcolo il prodotto scalare in L2([a,b])
					Ortog=0.;
					Ortog=0.5*(PSI[k][0]*PSI[j][0]+PSI[k][Nx-1]*PSI[j][Nx-1])*hx;
					for(i=1; i<Nx-1; i++)
					{
						Ortog+=PSI[k][i]*PSI[j][i]*hx;
					}
					cout<<"Prodotto scalare tra S"<<j<<" e S"<<k<<" : "<<Ortog<<endl;
				}
			}
			
		}
		
	}
	
	
		
} // fine del ciclo 0 degli stati eccitati

	cout<<"Fine programma."<<endl;	
 
 	return 0;
} 


// -----------------------------------------------------------------------------------------------------------------------------
// Definizioni di funzioni usate nel programma

// Funzione per il calcolo dell'energia
double Energia(double B_h[Nh],int n_stato, double PSI[N_stati][Nx], double x[Nx])
{
	int j=0;
	int i=0;
	int k=0;
	
    	double Energia=0.;
    	double f[Nx]={0.};
    	double Norm_f=0.;
    	double Coeff[n_stato]={0.}; 	// Per la ricerca degli stati eccitati
    	 
    
    	if (I_case==0)          // Caso quadratico
    	{
        	double B_hS[Nh]={0.};
        	B_hS[0]=1.;
        	B_hS[1]=2.;
        	B_hS[2]=3.;
        	B_hS[3]=4.;
        	Energia=0.;
        
       	for(j=0; j<Nh; j++)
        	{
           		Energia+=(B_h[j]-B_hS[j])*(B_h[j]-B_hS[j]);
        	}
        
    	} // Fine caso 0
    
    	if (I_case==1)          // Caso buca di potenziale con sinusoidi
    	{
        	 // Impongo tutte le componenti sinusoidali
         	for(j=0; j<Nh; j++)
        	{
           		for(i=0;  i<Nx; i++)
             		{
            			f[i]+=B_h[j]*sin((j+1)*x[i]*M_PI/x_max);
             		}
         	}
         	
         	//Normalizzo f
        	Norm_f=Norm(f);
         	for(i=0;i<Nx;i++)
         	{
            		 f[i]=f[i]/Norm_f;
        	}
    
    		//Stati eccitati
    		if(n_stato>0)
     		{
     		// Rendo f ortogonale agli autostati precedenti
     		for(k=0;k<n_stato;k++) 
     		{
		 	Coeff[k]=0.5*(f[0]*PSI[k][0]+f[Nx-1]*PSI[k][Nx-1])*hx;
		 	for(i=1;i<Nx-1;i++)
		 	{
		     		Coeff[k]+=f[i]*PSI[k][i]*hx;
			}
		}
		for(k=0;k<n_stato;k++)
     		{
			for(i=0; i<Nx; i++)
			{
		    		f[i]=f[i]-PSI[k][i]*Coeff[k];
			}
        	}

    	}
   
     	double F[Nx];
     	F[0]=0.;
     	F[Nx-1]=0.;
     	for(int i=1; i<Nx-1; i++)
     	{
     		F[i]=f[i]*(-0.5*(f[i+1]-2.*f[i]+f[i-1])+V(x[i])*f[i])/hx;
     	}	
    	//Integro
     	Energia=0.5*(F[0]+F[Nx-1]);
     	for(int i=1; i<Nx-1; i++)
     	{
        	Energia+=F[i];
     	}
        Norm_f=Norm(f);
        Energia=Energia/(Norm_f*Norm_f);
    } // end I_case==1

     return Energia;
 }
 

// funzione per il calcolo della derivata parziale
double DerivParziale(double B_h[Nh], int r, int n_stato, double PSI[N_stati][Nx],double x[Nx])
{
     double drp=0.;
     double B_h1[Nh]={0.};
     double B_h_1[Nh]={0.};
     double dB_h[Nh]={0.};
     dB_h[r]=h_der;
     
     for(int i=0; i<Nh; i++)
     {
         B_h1[i]=B_h[i]+dB_h[i];
         B_h_1[i]=B_h[i]-dB_h[i];
     }
     drp=(Energia(B_h1,n_stato,PSI,x)-Energia(B_h_1,n_stato,PSI,x))/(2.*h_der);
    return drp;
 }

// funzione potenziale V(x)
double V(double x)
{
    	double V=0.;
    	
        if((x<b)&&(x>a))
        {
            V=V0;
        }
        return V;
}

// funzione per il calcolo la norma in L2([x_min,x_max])
double Norm(double f[Nx])
{
    double Norm=0.;
    int i=0;
    Norm=0.5*(pow(f[0],2)+pow(f[Nx-1],2))*hx;
    for(i=1; i<Nx-1; i++)
    {
        Norm+=pow(f[i],2)*hx;
    }
    Norm=sqrt(Norm);
    return Norm;
}


// funzione per il calcolo il modulo di un vettore
double Modulo_v(double h[Nh])
{
    double Modulo_v=0.;
    int j=0;
    for(j=0; j<Nh; j++)
    {
        Modulo_v+=pow(h[j],2);
    }
    Modulo_v=sqrt(Modulo_v);
    return Modulo_v;
}
