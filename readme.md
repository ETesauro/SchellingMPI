# Schelling's Model

Progetto per l'esame di **Programmazione Concorrente e Parallela su Cloud**.

- Studente: **Tesauro Emmanuel**
- Matricola: **0522500988**
- MD5: **4bc390e61c1cd2a4001199d94b8d334d**
- Amazon EC2 instance types: **t2.2xlarge**

## Introduction

Il modello di segregazione di Schelling è un modello agent-based che dimostra che avere persone con una preferenza "lieve" all'interno del proprio gruppo nei confronti nel gruppo stesso, porta inevitabilmente ad una società segregata.\

L'obiettivo di tale progetto è quello di implementare il modello di segregazione di Schelling creando un programma scritto in linguaggio C ed utilizzando la libreria Open MPI.

Viene fornito un unico file (Schelling.c) che necessita di essere prima compilato e, poi, eseguito.\
Al suo interno è possibile definire la grandezza della matrice di partenza a proprio piacimento, in quanto l'algoritmo lavora con matrici di qualsiasi dimensione.

## Parallel Solution Description

<p style="color: #00aaff;"> DOING </p>
La soluzione sviluppata segue X passi ben definiti.\

- Il Master inizializza la matrice utilizzando le costanti definite all'interno del programma. Di default, questi valori sono:

  - numero di righe: **100**
  - numero di colonne: **100**
  - percentuale di agenti '<p style="color: #00aaff; display:inline">X</p>' all'interno della matrice: **30%**
  - percentuale di agenti '<p style="color: red; display:inline">O</p>' all'interno della matrice: **30%**
  - percentuale di soddisfazione degli agenti: **33.3%**

- Suddivisione della matrice per numero di righe tra i vari processi
- Esecuzione per un numero di volte pari MAX_STEP delle seguenti operazioni:

  - Scambio delle righe tra i processi adiacenti per il calcolo della soddisfazione (locale) degli agenti della sottomatrice (la prima riga con il processo precedente e l'ultima riga con il processo successivo).
  - Calcolo della soddisfazione (locale) di ogni agente della sottomatrice
  - Calcolo del numero di celle vuote di ogni sottomatrice e delle posizioni in cui sono situate
  - Scambio delle posizioni calcolate al punto precedente tra tutti i processi per l'assegnazione delle celle di destinazione in cui gli agenti possono spostarsi
  - Spostamento degli agenti

- Recupero della matrice finale per calcolare la soddisfazione globale

## Implementation details

<p style="color: orange;"> TODO </p>

## Execution instructions

<p style="color: orange;"> TODO </p>

## Correctness discussion

<p style="color: orange;"> TODO </p>

## Benchmarks

<p style="color: orange;"> TODO </p>
