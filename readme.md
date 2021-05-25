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

La soluzione sviluppata segue X passi ben definiti.

- Il Master inizializza la matrice utilizzando le costanti definite all'interno del programma. Di default, questi valori sono:

  - Numero di righe: **100**
  - Numero di colonne: **100**
  - Percentuale di agenti '**X**' all'interno della matrice: **30%**
  - Percentuale di agenti '**O**' all'interno della matrice: **30%**
  - Percentuale di soddisfazione degli agenti: **33.3%**

- Suddivisione della matrice per numero di righe tra i vari processi
- Esecuzione per un numero di volte pari MAX_STEP delle seguenti operazioni:

  - Scambio delle righe tra i processi adiacenti per il calcolo della soddisfazione (locale) degli agenti della sottomatrice (la prima riga con il processo precedente e l'ultima riga con il processo successivo).
  - Calcolo della soddisfazione (locale) di ogni agente della sottomatrice
  - Calcolo del numero di celle vuote di ogni sottomatrice e delle posizioni in cui sono situate
  - Scambio delle posizioni calcolate al punto precedente tra tutti i processi per l'assegnazione delle celle di destinazione in cui gli agenti possono spostarsi
  - Spostamento degli agenti

- Recupero della matrice finale per calcolare la soddisfazione globale

## Project structure

- src/
  - _Schelling_MPI.c_: file contenente il codice sorgente del programma
- files_out/
  - _Schelling_MPI.out_: file eseguibile del programma
  - _Schelling_MPI.html_: file generato al termine dell'esecuzione del programma che contiene la matrice finale.
- doc/
  - _mdb.min.css_: CSS utilizzato per la pagina HTML generata.

## Getting started

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

<p style="color: #00aaff;"> DOING </p>

## Execution instructions

<p style="color: orange;"> TODO </p>

## Correctness discussion

<p style="color: orange;"> TODO </p>

## Benchmarks

<p style="color: orange;"> TODO </p>
