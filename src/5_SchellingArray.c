#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#include "mpi.h"

#define TEST 0

// #region Matrice
#define ROWS 5     // Numero di righe della matrice
#define COLUMNS 5  // Numero di colonne della matrice
// #endregion

// #region Agenti
#define AGENT_X 'X'              // Agente X (BLUE)
#define AGENT_O 'O'              // Agente O (RED)
#define EMPTY ' '                // Casella vuota (' ')
#define AGENT_X_PERCENTAGE 40    // Percentuale di agenti X (blu) all'interno della matrice
#define AGENT_O_PERCENTAGE 40    // Percentuale di agenti O (rossi) all'interno della matrice
#define SATISFIED_PERCENTAGE 30  // Percentuale di soddisfazione di un agente
// #endregion

// #region Settings
#define ROOT 0
#define MAX_STEP 1
// #endregion

// #region Utils
#define RAND(min, max) ((rand() % (max)) + min)
#define PRINT_BLUE(str) printf("\033[1;34m%c\033[0m ", str);
#define PRINT_RED(str) printf("\033[1;31m%c\033[0m ", str);
// #endregion

typedef struct voidCell {
    int row_index;
    int column_index;
} voidCell;

// #region Function definitions

//DEBUG
void test_init_matrix(char matrix[ROWS * COLUMNS], int O_pct, int X_pct);
void test_print_matrix(int n1, int n2, int **a);
void print_int_matrix(int rows_size, int column_size, int *matrix);
void print_exchange_matrix(int rank, int world_size, int rows_size, int original_rows, int column_size, char *matrix);
//:DEBUG

// rs = rows_size = quante righe ha la matrice
// cs = columns_size = quante colonne ha la matrice
int init_matrix(char *mat, int, int);
int subdivide_matrix(int, int *, int *, int *);
void print_matrix(int rows, int columns, char *mat);
void exchange_rows(int, int, int rs, char *mat, MPI_Comm);
int *evaluate_move(int, int, int, int rs, char *mat);
int is_satisfied(int, int, int, int, int, int rs, char *mat);
void move(int, int, int, int rs, char *mat, int *);
int calculate_destination(int, int *, char *, char, int, int, int);
int try_move(int, int, int, int, char *, int, int, char *, int *);
/* 
void exchange_rows_ngh(int, int, int rs, int cs, char mat[rs][cs], MPI_Comm);
void calculate_satisfaction(char mat[ROWS][COLUMNS]);
void err_finish(int *, int *, int *); */
// #endregion

int main(int argc, char **argv) {
    int world_size, rank;
    double start_time, end_time;

    //Inizializzazione MPI
    MPI_Status status;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    //Esecuzione
    char matrix[ROWS * COLUMNS];  // Matrice di char ('X', 'O', ' ')
    int *displacements;           // Array che contiene i displacements per ogni processo (per la MPI_Scatterv)
    int *sendcounts;              // Array che contiene il numero di elementi (#righe_assegnate * #colonne) di un processo
    int *rows_per_process;        // Array che contiene il numero di righe assegnate ad ogni processo
    int *want_move = NULL;        // Array che indica quali agenti della sottomatrice vogliono muoversi

    start_time = MPI_Wtime();
    sendcounts = calloc(world_size, sizeof(int));
    displacements = calloc(world_size, sizeof(int));
    rows_per_process = calloc(world_size, sizeof(int));

    // Inizializzazione matrice
    if (rank == ROOT) {
        if (TEST)
            test_init_matrix(matrix, AGENT_O_PERCENTAGE, AGENT_X_PERCENTAGE);
        else
            init_matrix(matrix, AGENT_O_PERCENTAGE, AGENT_X_PERCENTAGE);
        print_matrix(ROWS, COLUMNS, matrix);
    }

    // Calcolo della porzione della matrice da assegnare a ciascun processo
    subdivide_matrix(world_size, displacements, sendcounts, rows_per_process);

    // Suddivisione delle righe tra i processi
    char sub_matrix[rows_per_process[rank] * COLUMNS];
    MPI_Scatterv(matrix, sendcounts, displacements, MPI_CHAR, sub_matrix, COLUMNS * ((ROWS / world_size) + 1), MPI_CHAR, ROOT, MPI_COMM_WORLD);

    // Calcolo di quante righe 'originali' ha il processo e di quante ne ha 'totali'
    int total_rows = rows_per_process[rank];
    int original_rows = total_rows - ((rank == 0 || rank == world_size - 1) ? 1 : 2);

    // Comincia l'esecuzione
    for (int i = 0; i < MAX_STEP; i++) {
        exchange_rows(rank, world_size, total_rows, sub_matrix, MPI_COMM_WORLD);             // Scambio le righe tra i processi
        want_move = evaluate_move(rank, world_size, original_rows, total_rows, sub_matrix);  // Vedo chi si vuole muovere

        sleep(rank);
        //printf("----------- %d -----------\n", rank);
        move(rank, world_size, original_rows, total_rows, sub_matrix, want_move);  // Lo muovo

        MPI_Barrier(MPI_COMM_WORLD);
    }

    /* int satisfied_agents = 0;
    for (int row = 0; row < original_rows; row++) {
        for (int column = 0; column < COLUMNS; column++) {
            if (sub_matrix[row * COLUMNS + column] != EMPTY)
                satisfied_agents = is_satisfied(rank, world_size, row, column, total_rows, original_rows, COLUMNS, sub_matrix) ? satisfied_agents + 1 : satisfied_agents;
        }
    }

    int total[world_size];
    MPI_Gather(&satisfied_agents, 1, MPI_INT, total, 1, MPI_INT, ROOT, MPI_COMM_WORLD);

    if (rank == ROOT) {
        int total_satisfied_agents = 0;
        for (int i = 0; i < world_size; i++) {
            total_satisfied_agents += total[i];
        }
        printf("\n\nIn tutto ci sono %d agenti soddisfatti\n", total_satisfied_agents);
        printf("Percentuale: %.f%%\n", (double)total_satisfied_agents / (ROWS * COLUMNS) * 100);
    } */

    end_time = MPI_Wtime();
    MPI_Finalize();

    if (rank == 0) {
        printf("\n🕒 Time in ms = %f\n", end_time - start_time);
        free(sendcounts);
        free(displacements);
        free(rows_per_process);
        free(want_move);
    }

    return 0;
}

// Funzione per inizializzare la matrice
int init_matrix(char *matrix, int O_pct, int X_pct) {
    int row, column, random;
    srand(time(NULL) + ROOT);

    if (ROWS <= 0 || COLUMNS <= 0) {
        printf("\033[1;31mERROR\033[0m! Matrix size must be lager than 0 * 0\n\n");
        return 0;
    }

    if (O_pct + X_pct >= 100) {
        printf("\033[1;31mERROR\033[0m! Red: %d%%, Blue: %d%%. The sum can't be >= 100.\n\n", O_pct, X_pct);
        return 0;
    }

    for (row = 0; row < ROWS; row++) {
        for (column = 0; column < COLUMNS; column++) {
            random = rand() % 100;

            if ((random >= 0) && (random < O_pct)) {
                *(matrix + (row * COLUMNS) + column) = AGENT_O;
            }

            if ((random >= O_pct) && (random < O_pct + X_pct)) {
                *(matrix + (row * COLUMNS) + column) = AGENT_X;
            }

            if ((random >= O_pct + X_pct) && (random < 100)) {
                *(matrix + (row * COLUMNS) + column) = EMPTY;
            }
        }
    }

    return 1;
}

// Funzione per calcolare la porzione di matrice per ogni processo
int subdivide_matrix(int world_size, int *displacements, int *sendcounts, int *rows_per_process) {
    if (world_size <= 0 || ROWS <= 0 || COLUMNS <= 0 || displacements == NULL || sendcounts == NULL || rows_per_process == NULL) {
        printf("\033[1;31mERROR\033[0m! Invalid input in subdivide_matrix.\n\n");
        return 0;
    }

    int divisione = (ROWS) / (world_size);
    int resto = (ROWS) % (world_size);
    int displacement = 0;

    for (int i = 0; i < world_size; i++) {
        sendcounts[i] = divisione * COLUMNS;
        rows_per_process[i] = divisione;
        if (resto > 0) {
            sendcounts[i] += COLUMNS;
            rows_per_process[i]++;
            resto--;
        }

        displacements[i] = displacement;
        displacement += sendcounts[i];

        // IDEA: Ho dato una riga in più al rank 0 e al rank (world_size - 1) e due righe in più ad ogni altro processo.
        rows_per_process[i] += (i == 0 || i == world_size - 1) ? 1 : 2;
    }

    return 1;
}

// Funzione per stampare la matrice
void print_matrix(int rows, int columns, char *mat) {
    int i;

    for (i = 0; i < rows * columns; i++) {
        printf("| ");
        if (*(mat + i) == AGENT_X) {
            PRINT_BLUE(*(mat + i));
        } else if (*(mat + i) == AGENT_O) {
            PRINT_RED(*(mat + i));
        } else
            printf("%c ", *(mat + i));

        if ((i + 1) % columns == 0 && i != 0)
            printf("|\n");
    }
    printf("\n");
}

// Funzione che scambia le righe di ogni processo con il proprio vicino
void exchange_rows(int rank, int world_size, int total_rows, char *sub_matrix, MPI_Comm communicator) {
    int neighbour_up, neighbour_down;
    MPI_Status status;
    MPI_Request request;

    neighbour_up = (rank + 1) % world_size;
    neighbour_down = (rank + world_size - 1) % world_size;

    int index = (rank == 0 || rank == world_size - 1) ? 1 : 2;

    if (rank != 0) {
        // Prima riga: a partire da 0 mando COLUMNS elementi
        MPI_Isend(sub_matrix, COLUMNS, MPI_CHAR, neighbour_down, 99, communicator, &request);
        MPI_Recv(sub_matrix + ((total_rows - index) * COLUMNS), COLUMNS, MPI_CHAR, neighbour_down, 99, communicator, &status);
    }

    if (rank != world_size - 1) {
        // Ultima riga: A partire da ((total_rows - index - 1) * COLUMNS) mando COLUMNS elementi
        MPI_Isend(sub_matrix + ((total_rows - index - 1) * COLUMNS), COLUMNS, MPI_CHAR, neighbour_up, 99, communicator, &request);
        MPI_Recv(sub_matrix + ((total_rows - index + ((rank == 0) ? 0 : 1)) * COLUMNS), COLUMNS, MPI_CHAR, neighbour_up, 99, communicator, &status);
    }
}

void send_rows(int rank, int world_size, int total_rows, char *sub_matrix, MPI_Comm communicator) {
    int neighbour_up, neighbour_down;
    MPI_Status status;

    neighbour_up = (rank + 1) % world_size;
    neighbour_down = (rank + world_size - 1) % world_size;

    int index = (rank == 0 || rank == world_size - 1) ? 1 : 2;

    if (rank != 0) {
        // Prima riga: a partire da 0 mando COLUMNS elementi
        MPI_Send(sub_matrix, COLUMNS, MPI_CHAR, neighbour_down, 99, communicator);
    }

    if (rank != world_size - 1) {
        // Ultima riga: A partire da ((total_rows - index - 1) * COLUMNS) mando COLUMNS elementi
        MPI_Send(sub_matrix + ((total_rows - index - 1) * COLUMNS), COLUMNS, MPI_CHAR, neighbour_up, 99, communicator);
    }
}

// Funzione per capire quale agente si deve spostare oppure no (in base alla soddisfazione)
int *evaluate_move(int rank, int world_size, int original_rows, int total_rows, char *sub_matrix) {
    // . 1: vuole spostarsi (non soddisfatto);
    // . 0: non si vuole spostare (soddisfatto);
    // . -1: l'elemento non esiste (libero per chi vuole spostarsi)
    int *want_move = malloc(total_rows * COLUMNS * sizeof(int));

    // Scorro prima le mie righe
    for (int row = 0; row < original_rows; row++) {
        for (int column = 0; column < COLUMNS; column++) {
            if (sub_matrix[(row * COLUMNS) + column] != EMPTY)
                want_move[(row * COLUMNS) + column] = is_satisfied(rank, world_size, row, column, total_rows, original_rows, sub_matrix) ? 0 : 1;
            else
                want_move[(row * COLUMNS) + column] = -1;
        }
    }

    // Scorro quelle dei vicini
    for (int row = original_rows; row < total_rows; row++) {
        for (int column = 0; column < COLUMNS; column++) {
            if (sub_matrix[(row * COLUMNS) + column] != EMPTY)
                want_move[(row * COLUMNS) + column] = 0;
            else
                want_move[(row * COLUMNS) + column] = -1;
        }
    }

    return want_move;
}

// Funzione che controlla se un agente è soddisfatto (1: soddisfatto; 0: non soddisfatto)
int is_satisfied(int rank, int world_size, int current_row, int current_column, int total_rows, int original_rows, char *sub_matrix) {
    char current_agent = *(sub_matrix + (current_row * COLUMNS) + current_column);  // Agente corrente
    int neighbours_count = 8;                                                       // Numero massimo di vicini
    int similar_agents = 0;                                                         // Agenti vicini uguali all'agente corrente

    // Indice del vicino sinistro (-1 se non esiste -> agente al bordo sinistro della matrice)
    int left_index = (COLUMNS + current_column - 1 != COLUMNS - 1) ? (COLUMNS + current_column - 1) % COLUMNS : -1;
    // Indice del vicino destro (-1 se non esiste -> agente al bordo destro della matrice)
    int right_index = ((COLUMNS + current_column + 1) % COLUMNS != 0) ? (COLUMNS + current_column + 1) % COLUMNS : -1;

    int ngh_precedent_row, ngh_next_row;
    char neighbours[8];

    // Calcolo la posizione delle righe precedenti e successive
    if (rank == 0) {
        ngh_precedent_row = -1;
        ngh_next_row = total_rows - 1;
    } else if (rank == world_size - 1) {
        ngh_precedent_row = total_rows - 1;
        ngh_next_row = -1;
    } else {
        ngh_precedent_row = total_rows - 1 - 1;
        ngh_next_row = total_rows - 1;
    }

    // Riga precedente
    if (current_row != 0) {
        if (left_index != -1)  // L'elemento a sinistra esiste
            neighbours[0] = *(sub_matrix + (current_row - 1) + left_index);
        else  // L'elemento a sinistra non esiste
            neighbours[0] = '\0';
    } else {
        if (left_index != -1) {
            neighbours[0] = rank == 0 ? '\0' : *(sub_matrix + (ngh_precedent_row * COLUMNS) + left_index);
        } else
            neighbours[0] = '\0';
    }

    if (current_row != 0) {
        neighbours[1] = *(sub_matrix + (current_row - 1) + current_column);
    } else {
        neighbours[1] = rank == 0 ? '\0' : *(sub_matrix + (ngh_precedent_row * COLUMNS) + current_column);
    }

    if (current_row != 0) {
        if (right_index != -1)  // L'elemento a destra esiste
            neighbours[2] = *(sub_matrix + (current_row - 1) + right_index);
        else  // L'elemento a destra non esiste
            neighbours[2] = '\0';
    } else {
        if (right_index != -1)
            neighbours[2] = rank == 0 ? '\0' : *(sub_matrix + (ngh_precedent_row * COLUMNS) + right_index);
        else
            neighbours[2] = '\0';
    }

    // Riga corrente
    neighbours[3] = left_index != -1 ? *(sub_matrix + (current_row * COLUMNS) + left_index) : '\0';
    neighbours[4] = right_index != -1 ? *(sub_matrix + (current_row * COLUMNS) + right_index) : '\0';

    // Riga successiva
    if (current_row != original_rows - 1) {
        if (left_index != -1)  // L'elemento a sinistra esiste
            neighbours[5] = *(sub_matrix + ((current_row + 1) * COLUMNS) + left_index);
        else  // L'elemento a sinistra non esiste
            neighbours[5] = '\0';
    } else {
        if (left_index != -1) {
            //neighbours[5] = rank == 0 ? '\0' : sub_matrix[ngh_next_row][left_index];
            neighbours[5] = rank == world_size - 1 ? '\0' : *(sub_matrix + (ngh_next_row * COLUMNS) + left_index);
        } else
            neighbours[5] = '\0';
    }

    if (current_row != original_rows - 1) {
        neighbours[6] = *(sub_matrix + ((current_row + 1) * COLUMNS) + current_column);
    } else {
        neighbours[6] = rank == world_size - 1 ? '\0' : *(sub_matrix + (ngh_next_row * COLUMNS) + current_column);
    }

    if (current_row != original_rows - 1) {
        if (right_index != -1)  // L'elemento a destra esiste
            neighbours[7] = *(sub_matrix + ((current_row + 1) * COLUMNS) + right_index);
        else  // L'elemento a destra non esiste
            neighbours[7] = '\0';
    } else {
        if (right_index != -1)
            neighbours[7] = rank == world_size - 1 ? '\0' : *(sub_matrix + (ngh_next_row * COLUMNS) + right_index);
        else
            neighbours[7] = '\0';
    }

    for (int i = 0; i < 8; i++) {
        if (neighbours[i] == current_agent) {
            similar_agents++;
        } else if (neighbours[i] == '\0') {
            neighbours_count--;
        }
    }

    if ((((double)100 / neighbours_count) * similar_agents) >= SATISFIED_PERCENTAGE) {
        return 1;
    } else
        return 0;
}

// Funzione per capire se l'agente deve spostarsi
void move(int rank, int world_size, int original_rows, int total_rows, char *sub_matrix, int *want_move) {
    srand(time(NULL) + rank);

    // Scorro want_move per capire chi vuole muoversi
    for (int i = 0; i < original_rows * COLUMNS; i++)
        // Se l'agente vuole muoversi...
        if (want_move[i] == 1) {
            /* if (rank == 0)
                printf("sono [%d] NON mando - i: %d\n", rank, i); */
            char element = sub_matrix[i];
            int pos = i;

            int destination = calculate_destination(rank, want_move, sub_matrix, element, pos, original_rows, total_rows);

            if (destination != -1) {
                int row_ind = destination / COLUMNS;
                int col_ind = destination % COLUMNS;
                printf("Ho scelto di spostare \033[1;31m%c\033[0m da \x1b[36m[%d][%d]\x1b[0m a -> \x1b[32m[%d][%d]\x1b[0m\n", element, pos / COLUMNS, pos % COLUMNS, row_ind, col_ind);

                // - Creo la nuova riga (aggiornata)
                char new_row[COLUMNS];
                memcpy(new_row, sub_matrix + (row_ind * COLUMNS), COLUMNS);
                new_row[col_ind] = element;

                // - Debug: stampo
                printf("Nuova riga:\t| ");
                for (int i = 0; i < COLUMNS; i++) {
                    if (i == col_ind)
                        printf("\033[1;31m%c\033[0m | ", new_row[i]);
                    else
                        printf("%c | ", new_row[i]);
                }
                printf("\n");

                int done = try_move(rank, world_size, row_ind, pos, new_row, original_rows, total_rows, sub_matrix, want_move);
                /* if (done) {
                    sub_matrix[row_ind * COLUMNS + col_ind] = new_row[i];  // Sposta l'agente
                    sub_matrix[i] = EMPTY;                                 // Liberalo lo spazio nella sottomatrice

                    want_move[row_ind * COLUMNS + col_ind] = 0;  // Non rendere più disponibile lo spazio disponibile per altri
                    want_move[i] = -1;                           // Libera questo spazio precedente
                } */
            } else {
                //printf("Non ci sono caselle vuote per spostare \033[1;31m%c\033[0m da \x1b[36m[%d][%d]\x1b[0m\n", element, pos / COLUMNS, pos % COLUMNS);
            }
        }
        // Mando le mie righe esterne ai miei vicini perchè potrebbero servire a loro
        else {
            /* if (rank == 0)
                printf("sono [%d] e mando - i: %d\n", rank, i); */
            send_rows(rank, world_size, total_rows, sub_matrix, MPI_COMM_WORLD);
        }
}

int try_move(int rank, int world_size, int row_index, int source, char *new_row, int original_rows, int total_rows, char *sub_matrix, int *want_move) {
    char *rowToSendToLeft;
    char *rowToSendToRight;

    char riga_ricevuta_sinistra[COLUMNS];
    char riga_ricevuta_destra[COLUMNS];

    int offset = (rank == 0 || rank == world_size - 1) ? 1 : 2;
    int indice_riga_vicino_sinistro = -1;
    int indice_riga_vicino_destro = -1;

    int rank_vicino_sinistro = (rank + world_size - 1) % world_size;
    int rank_vicino_destro = (rank + 1) % world_size;

    MPI_Status status;
    MPI_Request request;

    // - Prendo la riga del il vicino precedente (rank != 0)
    if (rank != 0) {
        indice_riga_vicino_sinistro = total_rows - offset;
        rowToSendToLeft = sub_matrix + indice_riga_vicino_sinistro * COLUMNS;  // Mi prendo la riga del vicino sinistro
    }

    // - Prendo la riga del il vicino successivo (rank != world_size - 1)
    if (rank != world_size - 1) {
        indice_riga_vicino_destro = total_rows - offset + ((rank == 0) ? 0 : 1);
        rowToSendToRight = sub_matrix + indice_riga_vicino_destro * COLUMNS;  // Mi prendo la riga del vicino destro
    }

    if (row_index == indice_riga_vicino_sinistro) {
        printf("Devo spostare l'elemento nel vicino sinistro [%d]\n", rank_vicino_sinistro);
        rowToSendToLeft = new_row;
    } else if (row_index == indice_riga_vicino_destro) {
        printf("Devo spostare l'elemento nel vicino destro [%d]\n", rank_vicino_destro);
        rowToSendToRight = new_row;
    } else {
        printf("Devo spostare l'elemento nella mia sottomatrice\n");
    }

    // - Mando la riga del vicino destro (che sarà modificata oppure no, dipende dalla destinazione dell'agente)
    if (rank != 0) {
        MPI_Isend(rowToSendToLeft, COLUMNS, MPI_CHAR, rank_vicino_sinistro, 99, MPI_COMM_WORLD, &request);
        MPI_Recv(riga_ricevuta_sinistra, COLUMNS, MPI_CHAR, rank_vicino_sinistro, 99, MPI_COMM_WORLD, &status);
    }

    // - Mando le righe del vicino sinistro (che sarà modificata oppure no, dipende dalla destinazione dell'agente)
    if (rank != world_size - 1) {
        MPI_Isend(rowToSendToRight, COLUMNS, MPI_CHAR, rank_vicino_destro, 99, MPI_COMM_WORLD, &request);
        MPI_Recv(riga_ricevuta_destra, COLUMNS, MPI_CHAR, rank_vicino_destro, 99, MPI_COMM_WORLD, &status);
    }

    if (rank == 1) {
        print_matrix(1, COLUMNS, riga_ricevuta_sinistra);
        print_matrix(1, COLUMNS, riga_ricevuta_destra);
    }

    /* // - Comincio a fare i controlli
    int voglio_scrivere_nel_proprietario = -1;
    int voglio_scrivere_nel_vicino_sinistro = -1;
    int voglio_scrivere_nel_vicino_destro = -1;

    // - Controllo se ho scritto nella mia sottomatrice. Se sì, la aggiorno e rispondo ai vicini con -1
    if (row_index != indice_riga_vicino_destro && row_index != indice_riga_vicino_sinistro) {
        voglio_scrivere_nel_proprietario = 1;
        for (int i = 0; i < COLUMNS; i++) {
            if (sub_matrix[row_index * COLUMNS + i] != new_row[i]) {
                sub_matrix[row_index * COLUMNS + i] = new_row[i];  // Sposta l'agente
                sub_matrix[source] = EMPTY;                        // Liberalo lo spazio nella sottomatrice

                want_move[row_index * COLUMNS + i] = 0;  // Non rendere più disponibile lo spazio disponibile per altri
                want_move[source] = -1;                  // Libera questo spazio precedente
                break;
            }
        }
    }
    // Riga dei vicini -> Do la precedenza a chi scrive prima. Se tutti e due scrivono nella stessa cella, la precedenza ce l'ha il vicino sinistro
    else {
        for (int i = 0; i < COLUMNS; i++) {
            if (rank != 0) {
                if (riga_ricevuta_sinistra[i] != new_row[i]) {
                    voglio_scrivere_nel_vicino_sinistro = 1;
                    break;
                }
            }
            if (rank != world_size - 1) {
                if (riga_ricevuta_destra[i] != new_row[i]) {
                    voglio_scrivere_nel_vicino_destro = 1;
                    break;
                }
            }
        }
    } */

    /* int sin = -1, des = -1;
    if (rank != 0) {
        MPI_Isend(&voglio_scrivere_nel_vicino_sinistro, 1, MPI_INT, rank_vicino_sinistro, 99, MPI_COMM_WORLD, &request);
        MPI_Recv(&sin, 1, MPI_INT, rank_vicino_sinistro, 99, MPI_COMM_WORLD, &status);
    }

    if (rank != world_size - 1) {
        MPI_Isend(&voglio_scrivere_nel_vicino_destro, 1, MPI_INT, rank_vicino_destro, 99, MPI_COMM_WORLD, &request);
        MPI_Recv(&des, 1, MPI_INT, rank_vicino_destro, 99, MPI_COMM_WORLD, &status);
    }
    sleep(4);
    printf("[%d]\tsinistro vuole scrivere: %d\n\tdestro vuole scrivere: %d\n\n", rank, sin, des); */

    // - Gestione del conflitto
    return 1;
}

// Funzione per spostare l'agente
int calculate_destination(int rank, int *want_move, char *sub_matrix, char element, int pos, int original_rows, int total_rows) {
    int void_cells[total_rows * COLUMNS];
    int ind = 0;

    // Calcolo le celle vuote
    for (int i = 0; i < total_rows; i++)
        for (int j = 0; j < COLUMNS; j++)
            if (want_move[i * COLUMNS + j] == -1) {
                void_cells[ind] = i * COLUMNS + j;
                ind++;
            }

    // Se ci sono celle vuote...
    if (ind > 0) {
        // ...scelgo la cella di destinazione in modo casuale
        int rand_pos = RAND(0, ind);
        return void_cells[rand_pos];
    } else
        // Altrimenti non lo posso spostare
        return -1;
}

// #region : *************** DEBUG ***************
void test_init_matrix(char matrix[ROWS * COLUMNS], int O_pct, int X_pct) {
    *(matrix + 0) = AGENT_X;
    *(matrix + 1) = AGENT_O;
    *(matrix + 2) = AGENT_X;
    *(matrix + 3) = AGENT_X;
    *(matrix + 4) = AGENT_X;
    *(matrix + 5) = AGENT_X;
    *(matrix + 6) = AGENT_X;
    *(matrix + 7) = AGENT_X;
    *(matrix + 8) = EMPTY;
    *(matrix + 9) = AGENT_O;
    *(matrix + 10) = AGENT_X;
    *(matrix + 11) = AGENT_O;
    *(matrix + 12) = EMPTY;
    *(matrix + 13) = AGENT_O;
    *(matrix + 14) = AGENT_O;
    *(matrix + 15) = AGENT_O;
    *(matrix + 16) = EMPTY;
    *(matrix + 17) = AGENT_X;
    *(matrix + 18) = AGENT_O;
    *(matrix + 19) = AGENT_O;
    *(matrix + 20) = AGENT_X;
    *(matrix + 21) = AGENT_O;
    *(matrix + 22) = AGENT_O;
    *(matrix + 23) = AGENT_O;
    *(matrix + 24) = AGENT_O;
}

void print_int_matrix(int rows_size, int column_size, int *matrix) {
    int i, j;
    if (rows_size <= 0 || column_size <= 0 || matrix == NULL) {
        printf("\033[1;31mERROR\033[0m! Invalid input.\n");
        MPI_Abort(MPI_COMM_WORLD, 0);
        return;
    }
    for (i = 0; i < rows_size; i++) {
        for (j = 0; j < column_size; j++) {
            if (j == 0) printf("| ");
            printf("%d ", *(matrix + (i * column_size) + j));

            printf("| ");
        }
        printf("\n");
    }
}

void print_exchange_matrix(int rank, int world_size, int total_rows, int original_rows, int columns_size, char *matrix) {
    int i;
    if (total_rows <= 0 || columns_size <= 0 || matrix == NULL) {
        printf("\033[1;31mERROR\033[0m! Invalid input.\n");
        MPI_Abort(MPI_COMM_WORLD, 0);
        return;
    }

    int index = (rank == 0 || rank == world_size - 1) ? 1 : 2;

    printf("[%d]:\n", rank);
    if (rank != 0) {  // Riga precedente
        for (i = 0; i < columns_size; i++) {
            printf("| %c ", *(matrix + ((total_rows - (rank == world_size - 1 ? 1 : 2)) * columns_size + i)));
        }
        printf("|\n");
    }

    for (i = 0; i < original_rows * columns_size; i++) {
        printf("| ");
        if (*(matrix + i) == AGENT_X) {
            PRINT_BLUE(*(matrix + i));
        } else if (*(matrix + i) == AGENT_O) {
            PRINT_RED(*(matrix + i));
        } else
            printf("%c ", *(matrix + i));

        if ((i + 1) % columns_size == 0 && i != 0)
            printf("|\n");
    }

    if (rank != world_size - 1) {  // Riga successiva
        for (int i = 0; i < columns_size; i++) {
            printf("| %c ", *(matrix + ((total_rows - (rank == world_size - 1 ? 0 : 1)) * columns_size + i)));
        }
        printf("|\n");
    }
    printf("\n");
}
// #endregion : *************** DEBUG ***************