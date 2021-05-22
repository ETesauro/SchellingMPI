/**
 * * Schelling_1 V1
 * Gli agenti si possono spostare solo nelle celle vuote del processo corrente
*/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "mpi.h"

#define TEST 0

// #region Matrice
#define ROWS 10     // Numero di righe della matrice
#define COLUMNS 10  // Numero di colonne della matrice
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
#define MAX_STEP 100
// #endregion

// #region Utils
#define RAND(min, max) ((rand() % (max)) + min)
#define PRINT_BLUE(str) printf("\033[1;34m%c\033[0m ", str);
#define PRINT_RED(str) printf("\033[1;31m%c\033[0m ", str);

#define BLUE(string) "\033[1;34m" string "\x1b[0m"
#define RED(string) "\033[1;31m" string "\x1b[0m"
#define YELLOW(string) "\033[1;33m" string "\x1b[0m"
// #endregion

// #region Structs
typedef struct voidCell {
    int row_index;
    int column_index;
} voidCell;
// #endregion

// #region : *************** DEBUG ***************
void test_init_matrix(char matrix[ROWS][COLUMNS], int O_pct, int X_pct) {
    int row, column, random;

    *(*(matrix + 0) + 0) = AGENT_X;
    *(*(matrix + 0) + 1) = AGENT_O;
    *(*(matrix + 0) + 2) = AGENT_X;
    *(*(matrix + 0) + 3) = AGENT_X;
    *(*(matrix + 0) + 4) = AGENT_X;
    *(*(matrix + 1) + 0) = AGENT_X;
    *(*(matrix + 1) + 1) = AGENT_X;
    *(*(matrix + 1) + 2) = AGENT_X;
    *(*(matrix + 1) + 3) = AGENT_O;
    *(*(matrix + 1) + 4) = AGENT_O;
    *(*(matrix + 2) + 0) = AGENT_X;
    *(*(matrix + 2) + 1) = AGENT_O;
    *(*(matrix + 2) + 2) = EMPTY;
    *(*(matrix + 2) + 3) = AGENT_O;
    *(*(matrix + 2) + 4) = AGENT_O;
    *(*(matrix + 3) + 0) = AGENT_O;
    *(*(matrix + 3) + 1) = EMPTY;
    *(*(matrix + 3) + 2) = EMPTY;
    *(*(matrix + 3) + 3) = AGENT_O;
    *(*(matrix + 3) + 4) = AGENT_O;
    *(*(matrix + 4) + 0) = EMPTY;
    *(*(matrix + 4) + 1) = AGENT_O;
    *(*(matrix + 4) + 2) = AGENT_O;
    *(*(matrix + 4) + 3) = AGENT_X;
    *(*(matrix + 4) + 4) = AGENT_O;
}

void test_print_matrix(int n1, int n2, int **a) {
    //Haoli's code
    // TODO: Your code here. Don't forget to:
    int i, j;
    if (n1 <= 0 || n2 <= 0 || a == NULL) {
        printf("input unvalid");
        return;
    }
    for (i = 0; i < n1; i++) {
        printf("| ");
        for (j = 0; j < n2; j++)
            printf("%d | ", *(*(a + i) + j));
        printf("\n");
    }
}

void print_int_matrix(int rows_size, int column_size, int matrix[rows_size][column_size]) {
    int i, j;
    if (rows_size <= 0 || column_size <= 0 || matrix == NULL) {
        printf("\033[1;31mERROR\033[0m! Invalid input.\n");
        MPI_Abort(MPI_COMM_WORLD, 0);
        return;
    }
    for (i = 0; i < rows_size; i++) {
        for (j = 0; j < column_size; j++) {
            if (j == 0) printf("| ");
            printf("%d ", *(*(matrix + i) + j));

            printf("| ");
        }
        printf("\n");
    }
}

// #endregion: *************** DEBUG ***************

// Funzione per inizializzare la matrice
int init_matrix(char matrix[ROWS][COLUMNS], int O_pct, int X_pct) {
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
                *(*(matrix + row) + column) = AGENT_O;
            }

            if ((random >= O_pct) && (random < O_pct + X_pct)) {
                *(*(matrix + row) + column) = AGENT_X;
            }

            if ((random >= O_pct + X_pct) && (random < 100)) {
                *(*(matrix + row) + column) = EMPTY;
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

//Funzione per stampare la matrice
void print_matrix(int rows_size, int column_size, char matrix[rows_size][column_size]) {
    int i, j;
    if (rows_size <= 0 || column_size <= 0 || matrix == NULL) {
        printf("\033[1;31mERROR\033[0m! Invalid input.\n");
        MPI_Abort(MPI_COMM_WORLD, 0);
        return;
    }
    for (i = 0; i < rows_size; i++) {
        for (j = 0; j < column_size; j++) {
            if (j == 0) printf("| ");
            if (*(*(matrix + i) + j) == AGENT_X) {
                printf(BLUE("%c") " ", *(*(matrix + i) + j));
            } else if (*(*(matrix + i) + j) == AGENT_O) {
                printf(RED("%c") " ", *(*(matrix + i) + j));
            } else
                printf("%c ", *(*(matrix + i) + j));

            printf("| ");
        }
        printf("\n");
    }
}

// Funzione che scambia le righe di ogni processo con il proprio vicino
void exchange_rows(int rank, int world_size, int rows, int columns, char sub_matrix[rows][columns], MPI_Comm communicator) {
    int neighbour_up, neighbour_down;
    //int received_row[columns];
    MPI_Status status;
    MPI_Request request;

    neighbour_up = (rank + 1) % world_size;
    neighbour_down = (rank + world_size - 1) % world_size;

    int index = (rank == 0 || rank == world_size - 1) ? 1 : 2;

    char *first_row = sub_matrix[0];
    char *last_row = sub_matrix[rows - index - 1];

    if (rank != 0) {
        MPI_Isend(first_row, columns, MPI_CHAR, neighbour_down, 99, communicator, &request);
        MPI_Recv(sub_matrix[rows - index], columns, MPI_CHAR, neighbour_down, 99, communicator, &status);
    }

    if (rank != world_size - 1) {
        MPI_Isend(last_row, columns, MPI_CHAR, neighbour_up, 99, communicator, &request);
        MPI_Recv(sub_matrix[rows - index + ((rank == 0) ? 0 : 1)], columns, MPI_CHAR, neighbour_up, 99, communicator, &status);
    }
}

// Funzione che controlla se un agente è soddisfatto (1: soddisfatto; 0: non soddisfatto)
int is_satisfied(int rank, int world_size, int rows_size, int total_rows, int columns_size, int row, int column, char sub_matrix[rows_size][columns_size]) {
    int left_index, right_index;
    int ngh_precedent_row, ngh_next_row;
    char neighbours[8];

    char current_element = sub_matrix[row][column];
    int neighbours_count = 8;
    int similar = 0;

    left_index = columns_size + column - 1 != columns_size - 1 ? (columns_size + column - 1) % columns_size : -1;
    right_index = ((columns_size + column + 1) % columns_size != 0) ? (columns_size + column + 1) % columns_size : -1;

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
    if (row != 0) {
        if (left_index != -1)  // L'elemento a sinistra esiste
            neighbours[0] = sub_matrix[row - 1][left_index];
        else  // L'elemento a sinistra non esiste
            neighbours[0] = '\0';
    } else {
        if (left_index != -1) {
            neighbours[0] = rank == 0 ? '\0' : sub_matrix[ngh_precedent_row][left_index];
        } else
            neighbours[0] = '\0';
    }

    if (row != 0) {
        neighbours[1] = sub_matrix[row - 1][column];
    } else {
        neighbours[1] = rank == 0 ? '\0' : sub_matrix[ngh_precedent_row][column];
    }

    if (row != 0) {
        if (right_index != -1)  // L'elemento a destra esiste
            neighbours[2] = sub_matrix[row - 1][right_index];
        else  // L'elemento a destra non esiste
            neighbours[2] = '\0';
    } else {
        if (right_index != -1)
            neighbours[2] = rank == 0 ? '\0' : sub_matrix[ngh_precedent_row][right_index];
        else
            neighbours[2] = '\0';
    }

    // Riga corrente
    neighbours[3] = left_index != -1 ? sub_matrix[row][left_index] : '\0';
    neighbours[4] = right_index != -1 ? sub_matrix[row][right_index] : '\0';

    // Riga successiva
    if (row != rows_size - 1) {
        if (left_index != -1)  // L'elemento a sinistra esiste
            neighbours[5] = sub_matrix[row + 1][left_index];
        else  // L'elemento a sinistra non esiste
            neighbours[5] = '\0';
    } else {
        if (left_index != -1) {
            neighbours[5] = rank == world_size - 1 ? '\0' : sub_matrix[ngh_next_row][left_index];
        } else
            neighbours[5] = '\0';
    }

    if (row != rows_size - 1) {
        neighbours[6] = sub_matrix[row + 1][column];
    } else {
        neighbours[6] = rank == world_size - 1 ? '\0' : sub_matrix[ngh_next_row][column];
    }

    if (row != rows_size - 1) {
        if (right_index != -1)  // L'elemento a destra esiste
            neighbours[7] = sub_matrix[row + 1][right_index];
        else  // L'elemento a destra non esiste
            neighbours[7] = '\0';
    } else {
        if (right_index != -1)
            neighbours[7] = rank == world_size - 1 ? '\0' : sub_matrix[ngh_next_row][right_index];
        else
            neighbours[7] = '\0';
    }

    for (int i = 0; i < 8; i++) {
        if (neighbours[i] == current_element) {
            similar++;
        } else if (neighbours[i] == '\0')
            neighbours_count--;
    }

    if ((((double)100 / neighbours_count) * similar) >= SATISFIED_PERCENTAGE)
        return 1;
    else
        return 0;
}

// Funzione per capire quale agente si deve spostare oppure no
int **evaluate_move(int rank, int world_size, int rows_size, int columns_size, char sub_matrix[rows_size][columns_size], int total_rows) {
    // 1: vuole spostarsi (non soddisfatto);
    // 0: non si vuole spostare (soddisfatto);
    // -1: l'elemento non esiste (libero per chi vuole spostarsi)
    int **mat = (int **)malloc(sizeof(int *) * rows_size);
    for (int i = 0; i < rows_size; i++)
        mat[i] = (int *)malloc(sizeof(int) * columns_size);

    for (int row = 0; row < rows_size; row++) {
        for (int column = 0; column < columns_size; column++) {
            if (sub_matrix[row][column] != EMPTY)
                *(*(mat + row) + column) = is_satisfied(rank, world_size, rows_size, total_rows, columns_size, row, column, sub_matrix) ? 0 : 1;
            else
                *(*(mat + row) + column) = -1;
        }
    }

    return mat;
}

// Funzione per spostare l'agente
void move(int row, int column, int original_rows, int columns, char mat[original_rows][columns], int **want_move) {
    /* for (int i = 0; i < original_rows; i++) {
        for (int j = 0; j < columns; j++) {
            if (*(*(want_move + i) + j) == -1) {
                mat[i][j] = mat[row][column];  // Sposta l'agente
                mat[row][column] = EMPTY;      // Liberalo lo spazio nella sottomatrice

                *(*(want_move + i) + j) = 0;          // Non rendere più disponibile lo spazio disponibile per altri
                *(*(want_move + row) + column) = -1;  // Libera questo spazio precedente
                break;
            }
        }
    } */

    voidCell tmp[original_rows * columns];
    int ind = 0;
    for (int i = 0; i < original_rows; i++)
        for (int j = 0; j < columns; j++)
            if (*(*(want_move + i) + j) == -1) {
                voidCell temp = {i, j};
                tmp[ind] = temp;
                ind++;
            }

    if (ind > 0) {
        int rand_pos = RAND(0, ind);
        voidCell var = tmp[rand_pos];
        mat[var.row_index][var.column_index] = mat[row][column];  // Sposta l'agente
        mat[row][column] = EMPTY;                                 // Liberalo lo spazio nella sottomatrice

        *(*(want_move + var.row_index) + var.column_index) = 0;  // Non rendere più disponibile lo spazio disponibile per altri
        *(*(want_move + row) + column) = -1;                     // Libera questo spazio precedente
    }
}

// Funzione per capire se l'agente deve spostarsi
void check_to_move(int rank, int original_rows, int total_rows, int columns, char sub_matrix[total_rows][columns], int **want_move) {
    for (int row = 0; row < original_rows; row++)
        for (int column = 0; column < columns; column++) {
            // Sposta l'agente in una cella libera
            if (*(*(want_move + row) + column) == 1)
                move(row, column, original_rows, columns, sub_matrix, want_move);
        }
}

void calculate_total_satisfaction(int rank, int world_size, char matrix[ROWS][COLUMNS]) {
    int total_agents = 0;
    int satisfied_agents = 0;

    for (int i = 0; i < ROWS; i++) {
        for (int j = 0; j < COLUMNS; j++) {
            if (matrix[i][j] != EMPTY) {
                total_agents++;
                if (is_satisfied(rank, world_size, ROWS, ROWS, COLUMNS, i, j, matrix)) {
                    satisfied_agents++;
                }
            }
        }
    }

    float average = ((double)satisfied_agents / (double)total_agents) * 100;
    printf("\nAgenti totali: %d\n", total_agents);
    printf("Agenti soddisfatti: %d\n", satisfied_agents);
    printf("\n🟢 Percentuale di soddisfazione: %.3f%%\n", average);
}

// Funzione per terminare l'esecuzione in caso di problemi
void err_finish(int *sendcounts, int *displacements, int *rows_per_process) {
    free(sendcounts);
    free(displacements);
    free(rows_per_process);
    MPI_Abort(MPI_COMM_WORLD, 0);
    exit(0);
}

int main(int argc, char **argv) {
    int world_size, rank;
    double start_time, end_time;

    //Inizializzazione MPI
    MPI_Status status;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    //Esecuzione
    char matrix[ROWS][COLUMNS];  // Matrice di char ('X', 'O', ' ')
    int *displacements;          // Array che contiene i displacements per ogni processo (per la MPI_Scatterv)
    int *sendcounts;             // Array che contiene il numero di elementi (#righe_assegnate * #colonne) di un processo
    int *rows_per_process;       // Arrray che contiene il numero di righe assegnate ad ogni processo
    int **want_move = NULL;

    start_time = MPI_Wtime();
    sendcounts = calloc(world_size, sizeof(int));
    displacements = calloc(world_size, sizeof(int));
    rows_per_process = calloc(world_size, sizeof(int));
    assert(sendcounts != NULL && displacements != NULL && rows_per_process != NULL);

    if (rank == ROOT) {
        if (TEST)
            test_init_matrix(matrix, AGENT_O_PERCENTAGE, AGENT_X_PERCENTAGE);
        else if (!init_matrix(matrix, AGENT_O_PERCENTAGE, AGENT_X_PERCENTAGE))
            err_finish(sendcounts, displacements, rows_per_process);
        print_matrix(ROWS, COLUMNS, matrix);
    }

    if (!subdivide_matrix(world_size, displacements, sendcounts, rows_per_process))
        err_finish(sendcounts, displacements, rows_per_process);

    if (rank == ROOT) {
        printf("\n");
        for (int i = 0; i < world_size; i++) {
            printf("rpp[%d] = %d (%d + %d)\tsendcounts[%d] = %d\tdispls[%d] = %d\n", i, rows_per_process[i], rows_per_process[i] - ((i == 0 || i == world_size - 1) ? 1 : 2), (i == 0 || i == world_size - 1) ? 1 : 2, i, sendcounts[i], i, displacements[i]);
        }
        printf("Subdivide matrix: \x1b[32mDONE\033[0m\n\n");
    }

    char sub_matrix[rows_per_process[rank]][COLUMNS];
    MPI_Scatterv(matrix, sendcounts, displacements, MPI_CHAR, sub_matrix, COLUMNS * ((ROWS / world_size) + 1), MPI_CHAR, ROOT, MPI_COMM_WORLD);

    int total_rows = rows_per_process[rank];
    int original_rows = total_rows - ((rank == 0 || rank == world_size - 1) ? 1 : 2);
    for (int i = 0; i < MAX_STEP; i++) {
        exchange_rows(rank, world_size, total_rows, COLUMNS, sub_matrix, MPI_COMM_WORLD);
        want_move = evaluate_move(rank, world_size, original_rows, COLUMNS, sub_matrix, total_rows);
        check_to_move(rank, original_rows, total_rows, COLUMNS, sub_matrix, want_move);

        MPI_Barrier(MPI_COMM_WORLD);
    }

    MPI_Gatherv(sub_matrix, sendcounts[rank], MPI_CHAR, matrix, sendcounts, displacements, MPI_CHAR, ROOT, MPI_COMM_WORLD);

    end_time = MPI_Wtime();
    MPI_Finalize();

    // * Stampa matrice finale e calcolo della soddisfazione
    if (rank == ROOT) {
        printf("\n");
        print_matrix(ROWS, COLUMNS, matrix);
        calculate_total_satisfaction(rank, world_size, matrix);
        printf("\n🕒 Time in ms = %f\n", end_time - start_time);
    }
    free(sendcounts);
    free(displacements);
    free(rows_per_process);
    free(want_move);

    return 0;
}