/**
 * * Schelling_3 V1
 * Differisce dal programma 2_Schelling V1 per spostare gli agenti ovunque nella matrice e non solo localmente
*/

#include <assert.h>
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
#define AGENT_X 'X'                // Agente X (BLUE)
#define AGENT_O 'O'                // Agente O (RED)
#define EMPTY ' '                  // Casella vuota (' ')
#define AGENT_X_PERCENTAGE 40      // Percentuale di agenti X (blu) all'interno della matrice
#define AGENT_O_PERCENTAGE 40      // Percentuale di agenti O (rossi) all'interno della matrice
#define SATISFIED_PERCENTAGE 33.3  // Percentuale di soddisfazione di un agente
// #endregion

// #region Settings
#define ROOT 0
#define MAX_STEP 1
#define SEED 10
// #endregion

// #region Utils
#define RAND(min, max) ((rand() % (max)) + min)

#define BLUE(string) "\033[1;34m" string "\x1b[0m"
#define RED(string) "\033[1;31m" string "\x1b[0m"
#define YELLOW(string) "\033[1;33m" string "\x1b[0m"
// #endregion

// #region Structs
typedef struct voidCell {
    int row_index;
    int column_index;
} voidCell;

typedef struct moveAgent {
    int destination_row;
    int destination_column;
    char agent;
} moveAgent;
// #endregion

// #region Function Signatures
int init_matrix(char mat[ROWS][COLUMNS], int O_pct, int X_pct);                                                                                 // Funzione per inizializzare la matrice
int subdivide_matrix(int, int *, int *, int *);                                                                                                 // Funzione per suddividere la matrice tra i processi
void exchange_rows(int, int, int rs, char mat[rs][COLUMNS], MPI_Comm);                                                                          // Funzione per scambiare le righe di ogni processo con i propri vicini
int **evaluate_move(int, int, int or_rs, int tot_rs, char mat[or_rs][COLUMNS], int *);                                                          // Funzione per capire quale agente si deve spostare oppure no
int is_satisfied(int, int, int rs, int tot_rs, int, int, char sub_matrix[rs][COLUMNS]);                                                         // Funzione per controllare se un agente è soddisfatto (1: soddisfatto; 0: non soddisfatto)
voidCell *calculate_local_void_cells(int or_rs, char mat[or_rs][COLUMNS], int, int *);                                                          // Funzione per calcolare le mie celle vuote
voidCell *exchange_void_cells(int rank, int world_size, int number_of_local_void_cells, voidCell *local_void_cells, int *, MPI_Datatype, int);  // Funzione per unire tutte le celle vuote dei processi e restituire quelle di destinazione per il processo i-esimo
int calculate_source(int world_size, int *displacement, int *sendcounts, int row);                                                              // Funzione per calcolare a quale processo appartiene una determinata riga della matrice
void move(int, int, int rs, char mat[rs][COLUMNS], int **, voidCell *, int, int *, int *, MPI_Datatype);                                        // Funzione per spostare gli agenti
void synchronize(int, int ws, int *, int vc_ass, moveAgent **data, int or_rs, char mat[or_rs][COLUMNS], MPI_Datatype);                          // Funzione per sincronizzare gli spostamenti tra i processi
void calculate_total_satisfaction(int, int, char mat[ROWS][COLUMNS]);                                                                           // !: Funzione per calcolare la soddisfazione finale di tutti gli agenti della matrice

void define_voidCellType(MPI_Datatype *);   // Funzione per definire il tipo voidCell
void define_moveAgentType(MPI_Datatype *);  // Funzione per definire il tipo moveAgent

void print_matrix(int rs, int cs, char matrix[rs][cs]);         // Funzione per stampare la matrice
void save_to_file(int rs, int cs, char *fn, char mat[rs][cs]);  // Funzione per salvare l'output dell'ordinamento
void err_finish(int *, int *, int *);                           // Funzione per terminare l'esecuzione in caso di problemi

//DEBUG
void test_init_matrix(char matrix[ROWS][COLUMNS], int O_pct, int X_pct);
//:DEBUG
// #endregion

int main(int argc, char **argv) {
    int world_size, rank;
    double start_time, end_time;

    // - Inizializzazione MPI
    MPI_Status status;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // - Datatypes
    MPI_Datatype VOID_CELL_TYPE;
    define_voidCellType(&VOID_CELL_TYPE);

    MPI_Datatype MOVE_AGENT_TYPE;
    define_moveAgentType(&MOVE_AGENT_TYPE);

    // - Variabili
    char matrix[ROWS][COLUMNS];          // Matrice di char ('X', 'O', ' ')
    int *displacements = NULL;           // Array che contiene i displacements per ogni processo (per la MPI_Scatterv)
    int *sendcounts = NULL;              // Array che contiene il numero di elementi (#righe_assegnate * #colonne) di un processo
    int *rows_per_process = NULL;        // Array che contiene il numero di righe assegnate ad ogni processo
    int **want_move = NULL;              // Array che indica quali agenti della sottomatrice vogliono muoversi
    int unsatisfied_agents = 0;          // Numero di agenti insoddisfatti per ogni processo (per ogni iterazione)
    int number_of_local_void_cells = 0;  // Numero di celle vuote nella mia sottomatrice
    voidCell *local_void_cells = NULL;   // Array che contiene le celle vuote della mia sottomatrice
    int num_assigned_void_cells = 0;     // Numero di celle vuote che mi sono state assegnate
    voidCell *to_move;                   // Array che contiene le celle vuote che mi sono state assegnate dove poter spostare gli agenti

    start_time = MPI_Wtime();
    sendcounts = calloc(world_size, sizeof(int));
    displacements = calloc(world_size, sizeof(int));
    rows_per_process = calloc(world_size, sizeof(int));
    assert(sendcounts != NULL && displacements != NULL && rows_per_process != NULL);

    // * Inizializzazione matrice
    if (rank == ROOT) {
        if (TEST)
            test_init_matrix(matrix, AGENT_O_PERCENTAGE, AGENT_X_PERCENTAGE);
        else if (!init_matrix(matrix, AGENT_O_PERCENTAGE, AGENT_X_PERCENTAGE))
            err_finish(sendcounts, displacements, rows_per_process);
        print_matrix(ROWS, COLUMNS, matrix);
    }

    // * Calcolo della porzione della matrice da assegnare a ciascun processo
    if (!subdivide_matrix(world_size, displacements, sendcounts, rows_per_process))
        err_finish(sendcounts, displacements, rows_per_process);

    // * Suddivisione delle righe tra i processi
    char sub_matrix[rows_per_process[rank]][COLUMNS];
    MPI_Scatterv(matrix, sendcounts, displacements, MPI_CHAR, sub_matrix, COLUMNS * ((ROWS / world_size) + 1), MPI_CHAR, ROOT, MPI_COMM_WORLD);

    // * Calcolo di quante righe 'originali' ha il processo e di quante ne ha 'totali'
    int total_rows = rows_per_process[rank];
    int original_rows = total_rows - ((rank == 0 || rank == world_size - 1) ? 1 : 2);

    // * Comincia l'esecuzione
    for (int i = 0; i < MAX_STEP; i++) {
        exchange_rows(rank, world_size, total_rows, sub_matrix, MPI_COMM_WORLD);                                  // Scambio le righe tra i processi
        want_move = evaluate_move(rank, world_size, original_rows, total_rows, sub_matrix, &unsatisfied_agents);  // Vedo chi si vuole spostare

        local_void_cells = calculate_local_void_cells(original_rows, sub_matrix, displacements[rank], &number_of_local_void_cells);  // Calcolo le mie celle vuote

        printf("Sono il processo %d e stampo l'array: [ ", rank);
        for (int i = 0; i < number_of_local_void_cells; i++) {
            printf("%d %d ", local_void_cells[i].row_index, local_void_cells[i].column_index);
        }
        printf("]\n");

        to_move = exchange_void_cells(rank, world_size, number_of_local_void_cells, local_void_cells, &num_assigned_void_cells, VOID_CELL_TYPE, unsatisfied_agents);  // Metto tutte le celle vuote insieme e restituisco l'array con le mie celle di destinazione

        move(rank, world_size, original_rows, sub_matrix, want_move, to_move, num_assigned_void_cells, displacements, sendcounts, MOVE_AGENT_TYPE);  // Sposto gli agenti

        MPI_Barrier(MPI_COMM_WORLD);
    }

    // * Recupero la matrice finale
    MPI_Gatherv(sub_matrix, sendcounts[rank], MPI_CHAR, matrix, sendcounts, displacements, MPI_CHAR, ROOT, MPI_COMM_WORLD);

    end_time = MPI_Wtime();
    MPI_Type_free(&VOID_CELL_TYPE);
    MPI_Type_free(&MOVE_AGENT_TYPE);
    MPI_Finalize();

    // * Stampa matrice finale e calcolo della soddisfazione
    if (rank == ROOT) {
        printf("\n");
        //print_matrix(ROWS, COLUMNS, matrix);
        //calculate_total_satisfaction(rank, world_size, matrix);
        save_to_file(ROWS, COLUMNS, "final.txt", matrix);
        printf("\n🕒 Time in ms = %f\n", end_time - start_time);
    }

    free(sendcounts);
    free(displacements);
    free(rows_per_process);
    free(want_move);

    return 0;
}

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

void exchange_rows(int rank, int world_size, int total_rows, char sub_matrix[total_rows][COLUMNS], MPI_Comm communicator) {
    int neighbour_up, neighbour_down;
    MPI_Status status;
    MPI_Request request;

    neighbour_up = (rank + 1) % world_size;
    neighbour_down = (rank + world_size - 1) % world_size;

    int index = (rank == 0 || rank == world_size - 1) ? 1 : 2;

    char *first_row = sub_matrix[0];
    char *last_row = sub_matrix[total_rows - index - 1];

    if (rank != 0) {
        MPI_Isend(first_row, COLUMNS, MPI_CHAR, neighbour_down, 99, communicator, &request);
        MPI_Recv(sub_matrix[total_rows - index], COLUMNS, MPI_CHAR, neighbour_down, 99, communicator, &status);
    }

    if (rank != world_size - 1) {
        MPI_Isend(last_row, COLUMNS, MPI_CHAR, neighbour_up, 99, communicator, &request);
        MPI_Recv(sub_matrix[total_rows - index + ((rank == 0) ? 0 : 1)], COLUMNS, MPI_CHAR, neighbour_up, 99, communicator, &status);
    }
}

int **evaluate_move(int rank, int world_size, int original_rows, int total_rows, char sub_matrix[original_rows][COLUMNS], int *unsatisfied_agents) {
    // .1: vuole spostarsi (non soddisfatto);
    // .0: non si vuole spostare (soddisfatto);
    // .-1: l'elemento non esiste (libero per chi vuole spostarsi)
    int **mat = (int **)malloc(sizeof(int *) * original_rows);
    for (int i = 0; i < original_rows; i++)
        mat[i] = (int *)malloc(sizeof(int) * COLUMNS);

    *unsatisfied_agents = 0;

    for (int row = 0; row < original_rows; row++) {
        for (int column = 0; column < COLUMNS; column++) {
            if (sub_matrix[row][column] != EMPTY) {
                *(*(mat + row) + column) = is_satisfied(rank, world_size, original_rows, total_rows, row, column, sub_matrix) ? 0 : 1;
                if (*(*(mat + row) + column) == 1)
                    *unsatisfied_agents += 1;
            } else
                *(*(mat + row) + column) = -1;
        }
    }

    return mat;
}

int is_satisfied(int rank, int world_size, int rows_size, int total_rows, int row, int column, char sub_matrix[rows_size][COLUMNS]) {
    int left_index, right_index;
    int ngh_precedent_row, ngh_next_row;
    char neighbours[8];

    char current_element = sub_matrix[row][column];
    int neighbours_count = 8;
    int similar = 0;

    left_index = COLUMNS + column - 1 != COLUMNS - 1 ? (COLUMNS + column - 1) % COLUMNS : -1;
    right_index = ((COLUMNS + column + 1) % COLUMNS != 0) ? (COLUMNS + column + 1) % COLUMNS : -1;

    // * Calcolo la posizione delle righe precedenti e successive
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

    // * Riga precedente
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

    // * Riga corrente
    neighbours[3] = left_index != -1 ? sub_matrix[row][left_index] : '\0';
    neighbours[4] = right_index != -1 ? sub_matrix[row][right_index] : '\0';

    // * Riga successiva
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

voidCell *calculate_local_void_cells(int original_rows, char sub_matrix[original_rows][COLUMNS], int displacement, int *local_void_cells) {
    voidCell *void_cells = malloc(original_rows * COLUMNS * sizeof(voidCell));
    int ind = 0;

    // * Calcolo le celle vuote
    for (int i = 0; i < original_rows; i++)
        for (int j = 0; j < COLUMNS; j++)
            if (sub_matrix[i][j] == EMPTY) {
                voidCell temp = {(displacement / COLUMNS) + i, j};
                void_cells[ind] = temp;
                ind++;
            }

    void_cells = realloc(void_cells, ind * sizeof(voidCell));
    *local_void_cells = ind;

    return void_cells;
}

voidCell *exchange_void_cells(int rank, int world_size, int number_of_local_void_cells, voidCell *local_void_cells, int *number_of_void_cells_to_return, MPI_Datatype datatype, int unsatisfied_agents) {
    int number_of_global_void_cells[world_size];  // Array che contiene il numero di celle vuote per ogni processo
    int displacements[world_size];                // Displacements per la prima gather (displs[rank] == numero di celle vuote del rank)
    int number_of_total_void_cells = 0;           // Numero di celle vuote in TUTTA la matrice
    voidCell *global_void_cells;                  // Array con tutte le posizioni delle celle vuote - es. { [0][1], [3][2], [5][8] }
    int *void_cells_per_process;                  // Array che indica quante celle vuote vengono assegnate ad un processo
    int global_unsatisfied_agents[world_size];    // Array che contiene il numero degli agenti insoddisfatti per ogni processo

    global_void_cells = malloc(ROWS * COLUMNS * sizeof(voidCell));
    void_cells_per_process = malloc(ROWS * COLUMNS * sizeof(voidCell));

    // * Vedo prima gli altri quante celle vuote hanno
    MPI_Allgather(&number_of_local_void_cells, 1, MPI_INT, number_of_global_void_cells, 1, MPI_INT, MPI_COMM_WORLD);

    // * Calcolo i DISPLACEMENTS e quante celle vuote ho in tutta la matrice
    for (int i = 0; i < world_size; i++) {
        displacements[i] = i == 0 ? 0 : displacements[i - 1] + number_of_global_void_cells[i - 1];
        number_of_total_void_cells += number_of_global_void_cells[i];
    }

    // * Metto tutte le celle vuote della matrice insieme
    MPI_Allgatherv(local_void_cells, number_of_local_void_cells, datatype, global_void_cells, number_of_global_void_cells, displacements, datatype, MPI_COMM_WORLD);

    // * Metto insieme il numero di elementi insoddisfatti di ogni processo (per fare i calcoli)
    MPI_Allgather(&unsatisfied_agents, 1, MPI_INT, global_unsatisfied_agents, 1, MPI_INT, MPI_COMM_WORLD);

    // * Mischio l'array (uso lo STESSO seme per ogni processo in modo da avere lo shuffle uguale)
    srand(SEED);
    for (int i = 0; i < number_of_total_void_cells; i++) {
        int destination = RAND(0, number_of_total_void_cells);
        voidCell tmp = global_void_cells[destination];
        global_void_cells[destination] = global_void_cells[i];
        global_void_cells[i] = tmp;
    }

    // * Calcolo le posizioni vuote per me e le restituisco
    int resto = number_of_total_void_cells % world_size;
    int displacement = 0;
    for (int i = 0; i < world_size; i++) {
        int divisione = number_of_total_void_cells / world_size;
        void_cells_per_process[i] = divisione > global_unsatisfied_agents[i] ? global_unsatisfied_agents[i] : divisione;

        if (resto > 0 && !(divisione > global_unsatisfied_agents[i])) {
            void_cells_per_process[i]++;
            resto--;
        }

        displacements[i] = displacement;
        displacement += void_cells_per_process[i];
    }

    // * Restituisco al main
    *number_of_void_cells_to_return = void_cells_per_process[rank];  // Setto il numero di celle vuote che mi sono state assegnate
    voidCell *toReturn = malloc(sizeof(voidCell) * void_cells_per_process[rank]);
    MPI_Scatterv(global_void_cells, void_cells_per_process, displacements, datatype, toReturn, void_cells_per_process[rank], datatype, ROOT, MPI_COMM_WORLD);

    return toReturn;
}

int calculate_source(int world_size, int *displacement, int *sendcounts, int row) {
    int toReturn = 0;

    for (int tmp_rank = 0; tmp_rank < world_size; tmp_rank++) {
        int start = displacement[tmp_rank] / COLUMNS;
        int finish = (displacement[tmp_rank] / COLUMNS + sendcounts[tmp_rank] / COLUMNS);

        if (row >= start && row < finish) {
            toReturn = tmp_rank;
            break;
        }
    }
    return toReturn;
}

void move(int rank, int world_size, int original_rows, char sub_matrix[original_rows][COLUMNS], int **want_move, voidCell *destinations, int num_assigned_void_cells, int *displacements, int *sendcounts, MPI_Datatype move_agent_type) {
    int num_elems_to_send_to[world_size];                           // Array che contiene il numero di moveAgents da mandare al processo i-esimo
    memset(num_elems_to_send_to, 0, sizeof(num_elems_to_send_to));  // Per settare a 0 l'array che ho appena creato

    //moveAgent data[world_size][num_assigned_void_cells];  // Matrice che contiene sulle righe i processi e sulle colonne la cella di destinazione dell'agente che vuole spostarsi
    moveAgent **data;
    data = (moveAgent **)malloc(sizeof(moveAgent *) * world_size);
    for (int i = 0; i < world_size; i++)
        data[i] = (moveAgent *)malloc(sizeof(moveAgent) * num_assigned_void_cells);

    int used_void_cells_assigned = 0;  // Il numero delle celle vuote che mi sono state assegnate che ho usato. Non posso usarne più di quante me ne sono state assegnate

    for (int i = 0; i < original_rows && used_void_cells_assigned < num_assigned_void_cells; i++) {
        for (int j = 0; j < COLUMNS && used_void_cells_assigned < num_assigned_void_cells; j++) {
            // * Sposta l'agente in una cella libera
            if (want_move[i][j] == 1) {
                voidCell destination = destinations[used_void_cells_assigned];                                  // Prendo la destinazione
                int receiver = calculate_source(world_size, displacements, sendcounts, destination.row_index);  // Capisco chi è il destinatario

                if (destination.row_index == -1 || destination.column_index == -1)
                    return;

                // * Sono io -> lo sposto direttamente
                if (receiver == rank) {
                    int startRow = displacements[rank] / COLUMNS;    // Riga iniziale
                    int destRow = destination.row_index - startRow;  // Riga di destinazione

                    sub_matrix[destRow][destination.column_index] = sub_matrix[i][j];  // Sposta l'agente
                    sub_matrix[i][j] = EMPTY;                                          // Liberalo lo spazio nella sottomatrice

                    want_move[destRow][destination.column_index] = 0;  // Non rendere più disponibile lo spazio disponibile per altri
                    want_move[i][j] = -1;                              // Libera questo spazio precedente
                }
                // * Altrimenti impacchetto tutto
                else {
                    int startRow = displacements[receiver] / COLUMNS;  // Riga iniziale
                    int destRow = destination.row_index - startRow;    // Riga di destinazione

                    moveAgent var = {destRow, destination.column_index, sub_matrix[i][j]};
                    data[receiver][num_elems_to_send_to[receiver]] = var;
                    num_elems_to_send_to[receiver] += 1;

                    sub_matrix[i][j] = EMPTY;  // Liberalo lo spazio nella sottomatrice
                    want_move[i][j] = -1;      // Libera questo spazio precedente
                }

                // * Aggiorno il numero di celle vuote che ho usato
                used_void_cells_assigned++;
            }
        }
    }

    // * Sincronizzo tra tutti i processi
    synchronize(rank, world_size, num_elems_to_send_to, num_assigned_void_cells, data, original_rows, sub_matrix, move_agent_type);

    free(data);
}

void synchronize(int rank, int world_size, int *num_elems_to_send_to, int num_assigned_void_cells, moveAgent **data, int original_rows, char sub_matrix[original_rows][COLUMNS], MPI_Datatype move_agent_type) {
    int my_void_cell_used[world_size];
    MPI_Request requests1[world_size];
    MPI_Request requests2[world_size];

    moveAgent **moved_agents;  // Matrice degli agenti che ho ricevuto che devo aggiornare nella mia sottomatrice
    moved_agents = (moveAgent **)malloc(sizeof(moveAgent) * world_size);

    // * Dico ai num_elems_to_send_to quante celle loro ho usato
    for (int i = 0; i < world_size; i++) {
        if (i == rank) continue;

        MPI_Isend(&num_elems_to_send_to[i], 1, MPI_INT, i, 99, MPI_COMM_WORLD, &requests1[i]);  // Mando al processo i il numero di celle che ho usato
        MPI_Irecv(&my_void_cell_used[i], 1, MPI_INT, i, 99, MPI_COMM_WORLD, &requests1[i]);     // Ricevo dal processo i il numero di celle mie che lui ha usato
    }

    for (int i = 0; i < world_size; i++) {
        if (i == rank) continue;

        moveAgent toSend[num_elems_to_send_to[i]];

        // * Mando e ricevo al/dal processo i-esimo tutte le celle di destinazione dove devo scrivere/salvare i miei agenti
        for (int j = 0; j < num_elems_to_send_to[i]; j++) {
            moveAgent var = {
                data[i][j].destination_row,
                data[i][j].destination_column,
                data[i][j].agent};

            toSend[j] = var;

            //printf("var: [%d][%d - %c\n", var.destination_row, var.destination_column, var.agent);
        }

        MPI_Wait(&requests1[i], NULL);

        moved_agents[i] = (moveAgent *)malloc(sizeof(moveAgent) * my_void_cell_used[i]);
        MPI_Isend(toSend, num_elems_to_send_to[i], move_agent_type, i, 100, MPI_COMM_WORLD, &requests2[i]);
        MPI_Irecv(moved_agents[i], my_void_cell_used[i], move_agent_type, i, 100, MPI_COMM_WORLD, &requests2[i]);
    }

    for (int i = 0; i < world_size; i++) {
        if (i == rank) continue;

        MPI_Wait(&requests2[i], NULL);
    }

    /* for (int i = 0; i < world_size; i++) {
        if (i == rank) continue;

        MPI_Send(&num_elems_to_send_to[i], 1, MPI_INT, i, 99, MPI_COMM_WORLD);  // Mando al processo i il numero di celle che ho usato
        printf("[%d] 1 send- mando a %d\n", rank, i);
        MPI_Recv(&my_void_cell_used[i], 1, MPI_INT, i, 99, MPI_COMM_WORLD, NULL);  // Ricevo dal processo i il numero di celle mie che lui ha usato
        printf("[%d] 1 receive - ho ricevuto da %d\n", rank, i);

        moveAgent toSend[num_elems_to_send_to[i]];

        // * Mando e ricevo al/dal processo i-esimo tutte le celle di destinazione dove devo scrivere/salvare i miei agenti
        for (int j = 0; j < num_elems_to_send_to[i]; j++) {
            moveAgent var = {
                data[i][j].destination_row,
                data[i][j].destination_column,
                data[i][j].agent};

            toSend[j] = var;
        }

        moved_agents[i] = (moveAgent *)malloc(sizeof(moveAgent) * my_void_cell_used[i]);
        MPI_Send(toSend, num_elems_to_send_to[i], move_agent_type, i, 100, MPI_COMM_WORLD);
        printf("[%d] 2 send- mando a %d\n", rank, i);
        MPI_Recv(moved_agents[i], my_void_cell_used[i], move_agent_type, i, 100, MPI_COMM_WORLD, NULL);
        printf("[%d] 2 receive - ho ricevuto da %d\n", rank, i);
    } */

    // * Scrivo gli agenti 'nuovi' nelle celle di destinazione
    for (int i = 0; i < world_size; i++) {
        if (i == rank) continue;

        for (int k = 0; k < my_void_cell_used[i]; k++) {
            /* if (rank == 0)
                printf("dest: [%d][%d] - %c \t", moved_agents[i][k].destination_row, moved_agents[i][k].destination_column, moved_agents[i][k].agent); */
            /* moveAgent destination = moved_agents[i][k];
            sub_matrix[destination.destination_row][destination.destination_column] = destination.agent; */
        }
    }

    // * Dealloco la matrice
    for (int i = 0; i < world_size; i++) {
        if (i == rank) continue;
        free(moved_agents[i]);
    }
    free(moved_agents);
}

void define_voidCellType(MPI_Datatype *VOID_CELL_TYPE) {
    // - Void Cell Type - es. {int row: 3, int col: 4}
    int vc_nitems = 2;                // Numero di items
    int vc_block_length[2] = {1, 1};  // Lunghezza dei blocchi
    // Displacements
    MPI_Aint vc_offsets[2];
    voidCell cella_vuota;
    MPI_Aint vc_base_address;
    MPI_Get_address(&cella_vuota, &vc_base_address);
    MPI_Get_address(&cella_vuota.row_index, &vc_offsets[0]);
    MPI_Get_address(&cella_vuota.column_index, &vc_offsets[1]);
    vc_offsets[0] = MPI_Aint_diff(vc_offsets[0], vc_base_address);
    vc_offsets[1] = MPI_Aint_diff(vc_offsets[1], vc_base_address);
    //: Displacements
    MPI_Datatype vc_types[2] = {MPI_INT, MPI_INT};
    MPI_Type_create_struct(vc_nitems, vc_block_length, vc_offsets, vc_types, VOID_CELL_TYPE);
    MPI_Type_commit(VOID_CELL_TYPE);
}

void define_moveAgentType(MPI_Datatype *MOVE_AGENT_TYPE) {
    int ma_nitems = 3;                   // Numero di items
    int ma_block_length[3] = {1, 1, 1};  // Lunghezza dei blocchi
    // Displacements
    MPI_Aint ma_offsets[3];
    moveAgent move_agent;
    MPI_Aint ma_base_address;
    MPI_Get_address(&move_agent, &ma_base_address);
    MPI_Get_address(&move_agent.destination_row, &ma_offsets[0]);
    MPI_Get_address(&move_agent.destination_column, &ma_offsets[1]);
    MPI_Get_address(&move_agent.agent, &ma_offsets[2]);
    ma_offsets[0] = MPI_Aint_diff(ma_offsets[0], ma_base_address);
    ma_offsets[1] = MPI_Aint_diff(ma_offsets[1], ma_base_address);
    ma_offsets[2] = MPI_Aint_diff(ma_offsets[2], ma_base_address);
    //: Displacements
    MPI_Datatype ma_types[3] = {MPI_INT, MPI_INT, MPI_CHAR};
    MPI_Type_create_struct(ma_nitems, ma_block_length, ma_offsets, ma_types, MOVE_AGENT_TYPE);
    MPI_Type_commit(MOVE_AGENT_TYPE);
}

void calculate_total_satisfaction(int rank, int world_size, char matrix[ROWS][COLUMNS]) {
    int total_agents = 0;
    int satisfied_agents = 0;
    moveAgent not_satisfied_agents[ROWS * COLUMNS];
    int index = 0;

    for (int i = 0; i < ROWS; i++) {
        for (int j = 0; j < COLUMNS; j++) {
            if (matrix[i][j] != EMPTY) {
                total_agents++;
                if (is_satisfied(rank, world_size, ROWS, ROWS, i, j, matrix)) {
                    satisfied_agents++;
                } else {
                    moveAgent var = {i, j, matrix[i][j]};
                    not_satisfied_agents[index] = var;
                    index++;
                }
            }
        }
    }

    float average = ((double)satisfied_agents / (double)total_agents) * 100;
    printf("\nAgenti totali: %d\n", total_agents);
    printf("Agenti soddisfatti: %d\n", satisfied_agents);
    if (index != 0) {
        printf("Agenti non soddisfatti: %d\n", index);
        printf("| ");
        for (int i = 0; i < index; i++) {
            printf(YELLOW("%c") ": [%d][%d] | ", matrix[not_satisfied_agents[i].destination_row][not_satisfied_agents[i].destination_column], not_satisfied_agents[i].destination_row, not_satisfied_agents[i].destination_column);
        }
    }
    printf("\n\n🟢 Percentuale di soddisfazione: %.3f%%\n", average);
}

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

void save_to_file(int rows_size, int column_size, char *file_name, char matrix[rows_size][column_size]) {
    FILE *doc = fopen(file_name, "w");
    int i, j;

    for (i = 0; i < rows_size; i++) {
        for (j = 0; j < column_size; j++) {
            if (j == 0) {
                fprintf(doc, "|");
                fprintf(doc, " ");
            }

            fprintf(doc, "%c", matrix[i][j]);
            fprintf(doc, " ");

            fprintf(doc, "|");
            fprintf(doc, " ");
        }
        fprintf(doc, "\n");
    }
    fclose(doc);
}

void err_finish(int *sendcounts, int *displacements, int *rows_per_process) {
    free(sendcounts);
    free(displacements);
    free(rows_per_process);
    MPI_Abort(MPI_COMM_WORLD, 0);
    exit(0);
}

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
// #endregion: *************** DEBUG ***************