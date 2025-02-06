# README - Programma per il Calcolo delle Energie e degli Autostati in una Buca di Potenziale

## Introduzione Teorica
Questo programma calcola le energie e gli autostati di una particella in una buca di potenziale infinita con un gradino di potenziale interno. Il problema viene risolto utilizzando il metodo del gradiente coniugato per trovare i parametri ottimali delle autofunzioni.

Il sistema fisico considerato è una buca di potenziale delimitata da pareti impenetrabili, con una discontinuità rappresentata da un gradino di potenziale di altezza `V0`. La funzione d'onda della particella viene approssimata come una combinazione lineare di sinusoidi e ottimizzata per trovare gli autostati.

## Parametri di Inizializzazione
Il programma utilizza una serie di parametri definiti all'inizio del codice:

- **Configurazione del problema:**
  - `I_case = 1`: Seleziona il tipo di potenziale (0: quadratico, 1: buca di potenziale con sinusoidi).
  - `N_stati = 3`: Numero di autostati da calcolare.
  - `Nh = 15`: Numero di parametri per descrivere gli stati quantistici.
  - `V0 = 0.1`: Altezza del gradino di potenziale (in unità di Hartree).

- **Definizione dell'asse x:**
  - `x_min = 0.0`, `x_max = 10.0`: Estremi della buca di potenziale.
  - `a = 3.0`, `b = 7.0`: Posizione del gradino di potenziale.
  - `Nx = 200`: Numero di punti della discretizzazione dell'asse x.
  - `hx = (x_max - x_min) / (Nx - 1)`: Spaziatura tra i punti x.

- **Parametri per il metodo numerico:**
  - `h_der = 0.001`: Passo per la derivata numerica.
  - `toll1 = 1e-7`: Tolleranza per la ricerca del minimo lungo `u`.
  - `toll2 = 1e-8`: Tolleranza per il gradiente coniugato.
  - `LL = 0.01`: Incremento per il minimo lungo `u`.

## Struttura del Programma
Il programma è organizzato nei seguenti passi principali:

1. **Inizializzazione:**
   - Lettura dei nomi dei file di output.
   - Definizione dell'asse x e delle variabili necessarie.

2. **Ciclo sugli stati quantistici:**
   - Per ciascun autostato (`n_stato`):
     - Inizializzazione dei parametri.
     - Applicazione del metodo del gradiente coniugato per minimizzare l'energia.
     - Normalizzazione dell'autostato ottenuto.
     - Salvataggio dell'autostato e dell'energia corrispondente.
     - Verifica dell'ortogonalità con gli stati precedenti.

3. **Metodo del Gradiente Coniugato:**
   - Calcolo del gradiente dell'energia rispetto ai parametri.
   - Aggiornamento dei parametri lungo la direzione del gradiente.
   - Iterazione fino alla convergenza secondo le tolleranze definite (`toll1` e `toll2`).

4. **Calcolo delle Energie e delle Funzioni d'Onda:**
   - La funzione d'onda viene rappresentata come combinazione di sinusoidi.
   - L'energia viene calcolata come valore atteso dell'operatore hamiltoniano.
   - Gli stati eccitati vengono resi ortogonali ai precedenti.

5. **Output dei Risultati:**
   - Stampa dell'energia e dei parametri finali.
   - Scrittura degli autostati su file.
   - Verifica dell'ortogonalità degli autostati.

## Funzioni Utilizzate
Il programma utilizza diverse funzioni per calcoli specifici:

- `double Energia(double B_h[Nh], int n_stato, double PSI[N_stati][Nx], double x[Nx])`:
  - Calcola l'energia di un dato stato quantico.

- `double DerivParziale(double B_h[Nh], int r, int n_stato, double PSI[N_stati][Nx], double x[Nx])`:
  - Calcola la derivata dell'energia rispetto a un parametro.

- `double V(double x)`:
  - Definisce il potenziale della buca con gradino.

- `double Norm(double f[Nx])`:
  - Calcola la norma L2 di una funzione d'onda.

## Esecuzione del Programma
Per compilare ed eseguire il programma:

```sh
 g++ -o calcolo_autostati programma.cpp -lm
 ./calcolo_autostati
```

L'output include le energie degli autostati e i file contenenti le funzioni d'onda corrispondenti.

## Conclusione
Questo programma implementa un metodo numerico efficiente per determinare gli autostati di una particella in una buca di potenziale con un gradino. Utilizzando il metodo del gradiente coniugato, ottimizza le funzioni d'onda per minimizzare l'energia e garantisce l'ortogonalità degli stati calcolati.

