# Schelling's Model

Progetto per l'esame di **Programmazione Concorrente e Parallela su Cloud**.

- Studente: **Tesauro Emmanuel**
- Matricola: **0522500988**
- MD5: **4bc390e61c1cd2a4001199d94b8d334d**
- Amazon EC2 instance types: **t2.2xlarge**

## Introduction

Il modello di segregazione di Schelling è un modello agent-based che dimostra che avere persone con una preferenza "lieve" all'interno del proprio gruppo nei confronti nel gruppo stesso, porta inevitabilmente ad una società segregata.

L'obiettivo di tale progetto è quello di implementare il modello di segregazione di Schelling creando un programma scritto in linguaggio C ed utilizzando la libreria Open MPI.

Viene fornito un unico file (Schelling.c) che necessita di essere prima compilato e, poi, eseguito.\
Al suo interno è possibile definire la grandezza della matrice di partenza a proprio piacimento, in quanto l'algoritmo lavora con matrici di qualsiasi dimensione.

## Parallel Solution Description

La soluzione sviluppata segue 8 passi ben definiti.

1. Il Master inizializza la matrice utilizzando le costanti definite all'interno del programma. Di default, questi valori sono:

   - Numero di righe: **100**
   - Numero di colonne: **100**
   - Percentuale di agenti '**X**' all'interno della matrice: **30%**
   - Percentuale di agenti '**O**' all'interno della matrice: **30%**
   - Percentuale di soddisfazione degli agenti: **33.3%**

2. Suddivisione della matrice per numero di righe tra i vari processi
3. Scambio delle righe tra i processi adiacenti per il calcolo della soddisfazione (locale) degli agenti della sottomatrice (la prima riga con il processo precedente e l'ultima riga con il processo successivo).
4. Calcolo della soddisfazione di ogni agente della sottomatrice
5. Calcolo del numero di celle vuote di ogni sottomatrice e delle posizioni in cui sono situate
6. Scambio delle posizioni calcolate al punto precedente tra tutti i processi per l'assegnazione delle celle di destinazione in cui gli agenti possono spostarsi
7. Spostamento degli agenti
8. Recupero della matrice finale per calcolare la soddisfazione globale

> Nota: I punti **3-7** vengono ripetuti per un numero di volte pari **MAX_STEP**.

## Project structure

- src/
  - _**Schelling_MPI.c**_: file contenente il codice sorgente del programma
- files_out/
  - _**Schelling_MPI.out**_: file eseguibile del programma
  - _**Schelling_MPI.html**_: file generato al termine dell'esecuzione del programma che contiene la matrice finale
- doc/
  - _**mdb.min.css**_: CSS utilizzato per la pagina HTML generata

## Execution instructions

Dalla **root** del progetto:

1. Compilare il programma con

   ```sh
   mpicc src/Schelling_MPI.c -o files_out/Schelling_MPI.out
   ```

2. Eseguire il programma con

   ```sh
   mpirun --allow-run-as-root -np X files_out/Schelling_MPI.out
   ```

   dove 'X' è un numero intero.

## Implementation details

### Inizializzazione della matrice

La matrice viene inizializzata dal Master in base ai parametri definiti all'interno del codice sorgente del programma.\
L'agente da inserire in una cella **[i][j]** della matrice viene determinato calcolando un **numero casuale _num_** tra 0 e 99 e si inserisce:

- **X** ->&emsp;se 0 &le; num &lt; AGENT_X_PERCENTAGE
- **O** ->&emsp;se AGENT_X_PERCENTAGE &le; num &lt; AGENT_X_PERCENTAGE + AGENT_O_PERCENTAGE
- **EMPTY** (' ') ->&emsp;se AGENT_X_PERCENTAGE + AGENT_O_PERCENTAGE &le; num &lt; 100

Le uniche regole da rispettare sono:

- La somma della percentuale di probabilità di avere un agente 'X' o 'O' all'interno della matrice **non deve superare 99**
- Non è possibile creare una matrice con un numero di righe **maggiore** del numero di processi con cui si decide di eseguire il programma

### Suddivisione del carico di lavoro

Tutti i processi calcolano quante righe della matrice dovranno ricevere. Questo significa che ogni processo si occuperà di **'numero di righe ricevute &#215; numero di colonne'** elementi.

In base al numero di righe della matrice di partenza e al numero di processi con cui si esegue il programma, viene effettuata la divisione tra questi due valori e se il '**resto &gt; 0**' viene assegnata una riga in più al processo in questione.

Inoltre, sono state assegnate **ulteriori 2 righe** ai processi con **'rank &gt; 0'** e **'rank &lt; world_size-1'** e solo **'una ulteriore riga'** ai processi con **'rank == 0'** e **'rank == world_size - 1'** per 'ospitare' le righe dei processi adiacenti.

> Nota: I processi con **'rank == 0'** e **'rank == world_size - 1'** ricevono una singola riga in più perchè si è scelto di non implementare la matrice come un array circolare e, quindi, questi processi non sono vicini.

### Scambio delle righe tra processi adiacenti

<p style="color: orange;"> TODO </p>

### Calcolo della soddisfazione (locale)

<p style="color: orange;"> TODO </p>

### Spostamento degli agenti

<p style="color: orange;"> TODO </p>

## Correctness discussion

<p style="color: orange;"> TODO </p>

## Benchmarks

<p style="color: orange;"> TODO </p>
