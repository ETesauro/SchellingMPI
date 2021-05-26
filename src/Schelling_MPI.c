/**
 * * Schelling_3 V1
 * Differisce dal programma 2_Schelling V1 per spostare gli agenti ovunque nella matrice e non solo localmente
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
#define AGENT_X 'X'                // Agente X (BLUE)
#define AGENT_O 'O'                // Agente O (RED)
#define EMPTY ' '                  // Casella vuota (' ')
#define AGENT_X_PERCENTAGE 30      // Percentuale di agenti X (blu) all'interno della matrice
#define AGENT_O_PERCENTAGE 30      // Percentuale di agenti O (rossi) all'interno della matrice
#define SATISFIED_PERCENTAGE 33.3  // Percentuale di soddisfazione di un agente
// #endregion

// #region Settings
#define ROOT 0
#define MAX_STEP 100
#define SEED 15
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
int init_matrix(char *, int, int);                                                     // Funzione per inizializzare la matrice
int subdivide_matrix(int, int *, int *, int *);                                        // Funzione per suddividere la matrice tra i processi
void exchange_rows(int, int, int, char *, MPI_Comm);                                   // Funzione per scambiare le righe di ogni processo con i propri vicini
int *evaluate_move(int, int, int, int, char *, int *);                                 // Funzione per capire quale agente si deve spostare oppure no
int is_satisfied(int, int, int, int, int, int, char *);                                // Funzione per controllare se un agente è soddisfatto (1: soddisfatto; 0: non soddisfatto)
voidCell *calculate_local_void_cells(int, char *, int, int *);                         // Funzione per calcolare le mie celle vuote
voidCell *exchange_void_cells(int, int, int, voidCell *, int *, MPI_Datatype, int);    // Funzione per unire tutte le celle vuote dei processi e restituire quelle di destinazione per il processo i-esimo
int calculate_source(int, int *, int *, int);                                          // Funzione per calcolare a quale processo appartiene una determinata riga della matrice
void move(int, int, int, char *, int *, voidCell *, int, int *, int *, MPI_Datatype);  // Funzione per spostare gli agenti
void synchronize(int, int, int *, int, moveAgent **, int, char *, MPI_Datatype);       // Funzione per sincronizzare gli spostamenti tra i processi
void calculate_total_satisfaction(int, int, char *);                                   // Funzione per calcolare la soddisfazione finale di tutti gli agenti della matrice

void define_voidCellType(MPI_Datatype *);   // Funzione per definire il tipo voidCell
void define_moveAgentType(MPI_Datatype *);  // Funzione per definire il tipo moveAgent

void print_matrix(int, int, char *);          // Funzione per stampare la matrice
void save_to_file(int, int, char *, char *);  // Funzione per salvare l'output dell'ordinamento
void err_finish(int *, int *, int *);         // Funzione per terminare l'esecuzione in caso di problemi

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

    if (ROWS < world_size && rank == ROOT) {
        printf(RED("Error! You have more processes then ROWS. ") "ROWS: %d, world_size: %d\n", ROWS, world_size);
        MPI_Abort(MPI_COMM_WORLD, MPI_ERR_COUNT);
    }

    // - Datatypes
    MPI_Datatype VOID_CELL_TYPE;
    define_voidCellType(&VOID_CELL_TYPE);

    MPI_Datatype MOVE_AGENT_TYPE;
    define_moveAgentType(&MOVE_AGENT_TYPE);

    // - Variabili
    char *matrix = NULL;                  // Matrice di char ('X', 'O', ' ')
    char *sub_matrix = NULL;              // Sottomatrice che mi è stata assegnata
    int *displacements = NULL;            // Array che contiene i displacements per ogni processo (per la MPI_Scatterv)
    int *sendcounts = NULL;               // Array che contiene il numero di elementi (#righe_assegnate * #colonne) di un processo
    int *rows_per_process = NULL;         // Array che contiene il numero di righe assegnate ad ogni processo
    int *want_move = NULL;                // Array che indica quali agenti della sottomatrice vogliono muoversi
    int unsatisfied_agents = 0;           // Numero di agenti insoddisfatti per ogni processo (per ogni iterazione)
    int number_of_local_void_cells = 0;   // Numero di celle vuote nella mia sottomatrice
    voidCell *local_void_cells = NULL;    // Array che contiene le celle vuote della mia sottomatrice
    int number_of_destination_cells = 0;  // Numero di celle vuote che mi sono state assegnate
    voidCell *destinations = NULL;        // Array che contiene le celle vuote che mi sono state assegnate dove poter spostare gli agenti

    start_time = MPI_Wtime();
    sendcounts = calloc(world_size, sizeof(int));
    displacements = calloc(world_size, sizeof(int));
    rows_per_process = calloc(world_size, sizeof(int));
    assert(sendcounts != NULL && displacements != NULL && rows_per_process != NULL);

    // * Inizializzazione matrice
    if (rank == ROOT) {
        matrix = malloc(ROWS * COLUMNS * sizeof(char *));
        /* if (TEST)
            test_init_matrix(matrix, AGENT_O_PERCENTAGE, AGENT_X_PERCENTAGE);
        else  */
        if (!init_matrix(matrix, AGENT_O_PERCENTAGE, AGENT_X_PERCENTAGE))
            err_finish(sendcounts, displacements, rows_per_process);
        print_matrix(ROWS, COLUMNS, matrix);
    }

    // * Calcolo della porzione della matrice da assegnare a ciascun processo
    if (!subdivide_matrix(world_size, displacements, sendcounts, rows_per_process))
        err_finish(sendcounts, displacements, rows_per_process);

    // * Suddivisione delle righe tra i processi
    sub_matrix = malloc(rows_per_process[rank] * COLUMNS * sizeof(char *));
    MPI_Scatterv(matrix, sendcounts, displacements, MPI_CHAR, sub_matrix, rows_per_process[rank] * COLUMNS, MPI_CHAR, ROOT, MPI_COMM_WORLD);

    // * Calcolo di quante righe 'originali' ha il processo e di quante ne ha 'totali'
    int total_rows = rows_per_process[rank];
    int original_rows = total_rows - ((rank == 0 || rank == world_size - 1) ? 1 : 2);

    // * Comincia l'esecuzione
    for (int i = 0; i < MAX_STEP; i++) {
        exchange_rows(rank, world_size, original_rows, sub_matrix, MPI_COMM_WORLD);                               // Scambio le righe tra i processi
        want_move = evaluate_move(rank, world_size, original_rows, total_rows, sub_matrix, &unsatisfied_agents);  // Vedo chi si vuole spostare

        local_void_cells = calculate_local_void_cells(original_rows, sub_matrix, displacements[rank], &number_of_local_void_cells);                                            // Calcolo le mie celle vuote
        destinations = exchange_void_cells(rank, world_size, number_of_local_void_cells, local_void_cells, &number_of_destination_cells, VOID_CELL_TYPE, unsatisfied_agents);  // Metto tutte le celle vuote insieme e restituisco l'array con le mie celle di destinazione

        move(rank, world_size, original_rows, sub_matrix, want_move, destinations, number_of_destination_cells, displacements, sendcounts, MOVE_AGENT_TYPE);  // Sposto gli agenti
        MPI_Barrier(MPI_COMM_WORLD);
    }

    // * Recupero la matrice finale
    MPI_Gatherv(sub_matrix, sendcounts[rank], MPI_CHAR, matrix, sendcounts, displacements, MPI_CHAR, ROOT, MPI_COMM_WORLD);

    end_time = MPI_Wtime();
    MPI_Type_free(&VOID_CELL_TYPE);
    MPI_Type_free(&MOVE_AGENT_TYPE);
    MPI_Finalize();

    // * Stampa matrice finale e calcolo della soddisfazione totale
    if (rank == ROOT) {
        printf("\n");
        print_matrix(ROWS, COLUMNS, matrix);
        calculate_total_satisfaction(rank, world_size, matrix);
        //save_to_file(ROWS, COLUMNS, "../files_out/Schelling_MPI.html", matrix);
        printf("\n🕒 Time in ms = %f\n", end_time - start_time);
    }

    free(matrix);
    free(sendcounts);
    free(displacements);
    free(rows_per_process);
    free(want_move);

    return 0;
}

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

void exchange_rows(int rank, int world_size, int original_rows, char *sub_matrix, MPI_Comm communicator) {
    int neighbour_up, neighbour_down;
    MPI_Status status;
    MPI_Request request_up;
    MPI_Request request_down;

    neighbour_up = (rank + 1) % world_size;
    neighbour_down = (rank + world_size - 1) % world_size;

    int my_last_row_pos = (original_rows - 1) * COLUMNS;                           // Indice della mia ultima riga (quella che devo mandare al vicino superiore)
    int neighbour_down_row_pos = original_rows * COLUMNS;                          // Indice della riga dove andrò a 'salvare' l'ultima riga del vicino precedente
    int neighbour_up_row_pos = (original_rows + ((rank == 0) ? 0 : 1)) * COLUMNS;  //Indice della riga dove andrò a 'salvare' la prima riga del vicino superiore

    if (rank != 0) {
        MPI_Isend(sub_matrix, COLUMNS, MPI_CHAR, neighbour_down, 99, communicator, &request_up);
        MPI_Irecv(sub_matrix + neighbour_down_row_pos, COLUMNS, MPI_CHAR, neighbour_down, 99, communicator, &request_up);
    }

    if (rank != world_size - 1) {
        MPI_Isend(sub_matrix + my_last_row_pos, COLUMNS, MPI_CHAR, neighbour_up, 99, communicator, &request_down);
        MPI_Irecv(sub_matrix + neighbour_up_row_pos, COLUMNS, MPI_CHAR, neighbour_up, 99, communicator, &request_down);
    }

    if (rank != 0) {
        MPI_Wait(&request_up, NULL);
    }
    if (rank != world_size - 1) {
        MPI_Wait(&request_down, NULL);
    }
}

int *evaluate_move(int rank, int world_size, int original_rows, int total_rows, char *sub_matrix, int *unsatisfied_agents) {
    int *mat = (int *)malloc(original_rows * COLUMNS * sizeof(int));
    *unsatisfied_agents = 0;

    for (int i = 0; i < original_rows; i++) {
        for (int j = 0; j < COLUMNS; j++) {
            if (sub_matrix[i * COLUMNS + j] != EMPTY) {
                int satisfied = is_satisfied(rank, world_size, original_rows, total_rows, i * COLUMNS, j, sub_matrix);
                // is_satisfied = 1: vuole spostarsi (non soddisfatto)
                // is_satisfied = 0: non si vuole spostare (soddisfatto)
                mat[i * COLUMNS + j] = satisfied ? 0 : 1;

                // Se l'agente non è soddisfatto, incremento il numero dei miei agenti insoddisfatti (mi servirà dopo nel calcolo delle destinazioni possibili)
                if (mat[i * COLUMNS + j] == 1)
                    *unsatisfied_agents += 1;
            }
            // .-1: la cella è vuota (è libera per chi vuole spostarsi)
            else
                mat[i * COLUMNS + j] = -1;
        }
    }

    return mat;
}

int is_satisfied(int rank, int world_size, int rows_size, int total_rows, int row, int column, char *sub_matrix) {
    int left_index, right_index;
    int ngh_precedent_row, ngh_next_row;
    char neighbours[8];

    char current_element = sub_matrix[row + column];
    int neighbours_count = 8;
    int similar = 0;

    left_index = (COLUMNS + column - 1 != COLUMNS - 1) ? (COLUMNS + column - 1) % COLUMNS : -1;
    right_index = ((COLUMNS + column + 1) % COLUMNS != 0) ? (COLUMNS + column + 1) % COLUMNS : -1;

    // * Calcolo la posizione delle righe precedenti e successive
    if (rank == 0) {
        ngh_precedent_row = -1;
        ngh_next_row = total_rows * COLUMNS - COLUMNS;
    } else if (rank == world_size - 1) {
        ngh_precedent_row = total_rows * COLUMNS - COLUMNS;
        ngh_next_row = -1;
    } else {
        ngh_precedent_row = total_rows * COLUMNS - COLUMNS - COLUMNS;
        ngh_next_row = total_rows * COLUMNS - COLUMNS;
    }

    // * Riga precedente
    if (row != 0) {
        if (left_index != -1)  // L'elemento a sinistra esiste
            neighbours[0] = sub_matrix[row - COLUMNS + left_index];
        else  // L'elemento a sinistra non esiste
            neighbours[0] = '\0';
    } else {
        if (left_index != -1) {
            neighbours[0] = rank == 0 ? '\0' : sub_matrix[ngh_precedent_row + left_index];
        } else
            neighbours[0] = '\0';
    }

    if (row != 0) {
        neighbours[1] = sub_matrix[row - COLUMNS + column];
    } else {
        neighbours[1] = rank == 0 ? '\0' : sub_matrix[ngh_precedent_row + column];
    }

    if (row != 0) {
        if (right_index != -1)  // L'elemento a destra esiste
            neighbours[2] = sub_matrix[row - COLUMNS + right_index];
        else  // L'elemento a destra non esiste
            neighbours[2] = '\0';
    } else {
        if (right_index != -1)
            neighbours[2] = rank == 0 ? '\0' : sub_matrix[ngh_precedent_row + right_index];
        else
            neighbours[2] = '\0';
    }

    // * Riga corrente
    neighbours[3] = left_index != -1 ? sub_matrix[row + left_index] : '\0';
    neighbours[4] = right_index != -1 ? sub_matrix[row + right_index] : '\0';

    // * Riga successiva
    if (row != rows_size - 1) {
        if (left_index != -1)  // L'elemento a sinistra esiste
            neighbours[5] = sub_matrix[row + COLUMNS + left_index];
        else  // L'elemento a sinistra non esiste
            neighbours[5] = '\0';
    } else {
        if (left_index != -1) {
            neighbours[5] = rank == world_size - 1 ? '\0' : sub_matrix[ngh_next_row + left_index];
        } else
            neighbours[5] = '\0';
    }

    if (row != rows_size - 1) {
        neighbours[6] = sub_matrix[row + COLUMNS + column];
    } else {
        neighbours[6] = rank == world_size - 1 ? '\0' : sub_matrix[ngh_next_row + column];
    }

    if (row != rows_size - 1) {
        if (right_index != -1)  // L'elemento a destra esiste
            neighbours[7] = sub_matrix[row + COLUMNS + right_index];
        else  // L'elemento a destra non esiste
            neighbours[7] = '\0';
    } else {
        if (right_index != -1)
            neighbours[7] = rank == world_size - 1 ? '\0' : sub_matrix[ngh_next_row + right_index];
        else
            neighbours[7] = '\0';
    }

    for (int i = 0; i < 8; i++) {
        if (neighbours[i] == current_element)
            similar++;
        else if (neighbours[i] == '\0')
            neighbours_count--;
    }

    if ((((double)100 / neighbours_count) * similar) >= SATISFIED_PERCENTAGE)
        return 1;
    else
        return 0;
}

voidCell *calculate_local_void_cells(int original_rows, char *sub_matrix, int displacement, int *local_void_cells) {
    voidCell *void_cells;  // Array che contiene le mie celle vuote ([riga][colonna]) -> ([riga * COLUMNS + colonna])
    int ind = 0;           // Indice che indica quante celle vuote ho trovato. Alla fine lo assegno a 'local_void_cells'

    void_cells = malloc(original_rows * COLUMNS * sizeof(voidCell));

    // * Calcolo le celle vuote
    for (int i = 0; i < original_rows; i++)
        for (int j = 0; j < COLUMNS; j++)
            if (sub_matrix[i * COLUMNS + j] == EMPTY) {
                voidCell temp = {(displacement + i * COLUMNS), j};
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
    int divisione = number_of_total_void_cells / world_size;
    int resto = number_of_total_void_cells % world_size;
    int displacement = 0;
    for (int i = 0; i < world_size; i++) {
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
        int finish = (start + sendcounts[tmp_rank] / COLUMNS);

        if (row >= start && row < finish) {
            toReturn = tmp_rank;
            break;
        }
    }
    return toReturn;
}

void move(int rank, int world_size, int original_rows, char *sub_matrix, int *want_move, voidCell *destinations, int num_assigned_void_cells, int *displacements, int *sendcounts, MPI_Datatype move_agent_type) {
    int num_elems_to_send_to[world_size];  // Array che contiene il numero di moveAgents da mandare al processo i-esimo
    int used_void_cells_assigned = 0;      // Il numero delle celle vuote che mi sono state assegnate che ho usato. Non posso usarne più di quante me ne sono state assegnate
    moveAgent **data;                      // Matrice che contiene sulle righe i processi e sulle colonne la cella di destinazione dell'agente che vuole spostarsi

    memset(num_elems_to_send_to, 0, sizeof(num_elems_to_send_to));  // Setto a 0 gli elementi dell'array che ho appena creato
    data = (moveAgent **)malloc(sizeof(moveAgent *) * world_size);  // Alloco spazio per la matrice data
    for (int i = 0; i < world_size; i++)
        data[i] = (moveAgent *)malloc(sizeof(moveAgent) * num_assigned_void_cells);

    for (int i = 0; i < original_rows && used_void_cells_assigned < num_assigned_void_cells; i++) {
        for (int j = 0; j < COLUMNS && used_void_cells_assigned < num_assigned_void_cells; j++) {
            // * Sposta l'agente in una cella libera
            if (want_move[i * COLUMNS + j] == 1) {
                voidCell destination = destinations[used_void_cells_assigned];                                            // Prendo la destinazione
                int receiver = calculate_source(world_size, displacements, sendcounts, destination.row_index / COLUMNS);  // Capisco chi è il destinatario

                // * Sono io -> lo sposto direttamente
                if (receiver == rank) {
                    int startRow = displacements[rank];              // Mia riga iniziale
                    int destRow = destination.row_index - startRow;  // Mia riga di destinazione

                    sub_matrix[destRow + destination.column_index] = sub_matrix[i * COLUMNS + j];  // Sposta l'agente
                    sub_matrix[i * COLUMNS + j] = EMPTY;                                           // Liberalo lo spazio nella sottomatrice

                    want_move[destRow + destination.column_index] = 0;  // Non rendere più disponibile lo spazio disponibile per altri
                    want_move[i * COLUMNS + j] = -1;                    // Libera questo spazio precedente
                }
                // * Altrimenti impacchetto tutto
                else {
                    int startRow = displacements[receiver];          // Riga iniziale del destinatario
                    int destRow = destination.row_index - startRow;  // Riga di destinazione del destinatario

                    moveAgent var = {destRow, destination.column_index, sub_matrix[i * COLUMNS + j]};
                    data[receiver][num_elems_to_send_to[receiver]] = var;  // Setto, al processo 'receiver' la X-esima colonna con la cella di destinazione dell'agente
                    num_elems_to_send_to[receiver] += 1;                   //Aggiorno il numero di elementi che devo mandare al processo 'receiver'

                    sub_matrix[i * COLUMNS + j] = EMPTY;  // Liberalo lo spazio nella sottomatrice
                    want_move[i * COLUMNS + j] = -1;      // Libera questo spazio precedente
                }

                // * Aggiorno il numero di celle vuote che ho usato
                used_void_cells_assigned++;
            }
        }
    }

    // * Sincronizzo tutti i processi
    synchronize(rank, world_size, num_elems_to_send_to, num_assigned_void_cells, data, original_rows, sub_matrix, move_agent_type);

    free(data);
}

void synchronize(int rank, int world_size, int *num_elems_to_send_to, int num_assigned_void_cells, moveAgent **data, int original_rows, char *sub_matrix, MPI_Datatype move_agent_type) {
    int my_void_cell_used_by[world_size];  // Array che contiene in ogni cella il numero di elementi che il processo i-esimo vuole scrivere nelle mie celle della sottomatrice
    MPI_Request requests1[world_size];     // Array per le prime MPI_Irecv e MPI_Wait
    MPI_Request requests2[world_size];     // Array per le seconde MPI_Irecv e MPI_Wait
    moveAgent **moved_agents;              // Matrice degli agenti che ho ricevuto che devo aggiornare nella mia sottomatrice

    moved_agents = (moveAgent **)malloc(sizeof(moveAgent) * world_size);

    // * Dico ai num_elems_to_send_to quante celle loro ho usato
    for (int i = 0; i < world_size; i++) {
        if (i == rank) continue;

        MPI_Isend(&num_elems_to_send_to[i], 1, MPI_INT, i, 99, MPI_COMM_WORLD, &requests1[i]);  // Mando al processo i il numero di celle che ho usato
        MPI_Irecv(&my_void_cell_used_by[i], 1, MPI_INT, i, 99, MPI_COMM_WORLD, &requests1[i]);  // Ricevo dal processo i il numero di celle mie che lui ha usato
    }

    // * Mando e ricevo al/dal processo i-esimo tutte le celle di destinazione dove devo scrivere/salvare i miei agenti
    for (int i = 0; i < world_size; i++) {
        if (i == rank) continue;

        int number_of_elems_to_send;  // Numero di elementi da mandare al processo i-esimo
        moveAgent *elements_to_send;  // Array che contiene gli elementi da mandare al processo i-esimo

        number_of_elems_to_send = num_elems_to_send_to[i];
        elements_to_send = malloc(number_of_elems_to_send * sizeof(moveAgent));  // Alloco lo spazio per la i-esima riga

        for (int j = 0; j < number_of_elems_to_send; j++)
            elements_to_send[j] = data[i][j];

        MPI_Wait(&requests1[i], NULL);  // Aspetto che la prima Irecv riceva il numero di elementi che gli altri processi vogliono scrivere nelle mie celle della mia sottomatrice

        moved_agents[i] = (moveAgent *)malloc(my_void_cell_used_by[i] * sizeof(moveAgent));
        MPI_Isend(elements_to_send, number_of_elems_to_send, move_agent_type, i, 100, MPI_COMM_WORLD, &requests2[i]);
        MPI_Irecv(moved_agents[i], my_void_cell_used_by[i], move_agent_type, i, 100, MPI_COMM_WORLD, &requests2[i]);
    }

    // * Aspetto che le seconde MPI_Wait terminino
    for (int i = 0; i < world_size; i++) {
        if (i == rank) continue;

        MPI_Wait(&requests2[i], NULL);
    }

    // * Scrivo gli agenti 'nuovi' nelle celle di destinazione
    for (int i = 0; i < world_size; i++) {
        if (i == rank) continue;

        for (int k = 0; k < my_void_cell_used_by[i]; k++)
            sub_matrix[moved_agents[i][k].destination_row + moved_agents[i][k].destination_column] = moved_agents[i][k].agent;  // Scrivo l'agente nella mia cella vuota
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

void calculate_total_satisfaction(int rank, int world_size, char *matrix) {
    int total_agents = 0;
    int satisfied_agents = 0;
    moveAgent *not_satisfied_agents = malloc(ROWS * COLUMNS * sizeof(moveAgent));

    int index = 0;

    for (int i = 0; i < ROWS; i++)
        for (int j = 0; j < COLUMNS; j++)
            if (matrix[i * COLUMNS + j] != EMPTY) {
                total_agents++;
                if (is_satisfied(rank, world_size, ROWS, ROWS, i * COLUMNS, j, matrix)) {
                    satisfied_agents++;
                } else {
                    moveAgent var = {i * COLUMNS, j, matrix[i * COLUMNS + j]};
                    not_satisfied_agents[index] = var;
                    index++;
                }
            }

    float average = ((double)satisfied_agents / (double)total_agents) * 100;
    printf("\nAgenti totali: %d\n", total_agents);
    printf("Agenti soddisfatti: %d\n", satisfied_agents);
    if (index != 0) {
        printf("Agenti non soddisfatti: %d\n", index);
        /* printf("| ");
        for (int i = 0; i < index; i++)
            printf(YELLOW("%c") ": [%d][%d] | ", matrix[not_satisfied_agents[i].destination_row + not_satisfied_agents[i].destination_column], not_satisfied_agents[i].destination_row / COLUMNS, not_satisfied_agents[i].destination_column); */
    }
    printf("\n\n🟢 Percentuale di soddisfazione: %.3f%%\n", average);

    //free(not_satisfied_agents);
}

void print_matrix(int rows_size, int column_size, char *matrix) {
    int i;
    if (rows_size <= 0 || column_size <= 0 || matrix == NULL) {
        printf("\033[1;31mERROR\033[0m! Invalid input.\n");
        MPI_Abort(MPI_COMM_WORLD, 0);
        return;
    }

    for (i = 0; i < rows_size * column_size; i++) {
        printf("| ");

        if (*(matrix + i) == AGENT_X) {
            printf(BLUE("%c") " ", *(matrix + i));
        } else if (*(matrix + i) == AGENT_O) {
            printf(RED("%c") " ", *(matrix + i));
        } else
            printf("%c ", *(matrix + i));

        if ((i + 1) % column_size == 0 && i != 0)
            printf("|\n");
    }
}

void save_to_file(int rows_size, int column_size, char *file_name, char *matrix) {
    FILE *doc = fopen(file_name, "w");
    int i, j;

    fprintf(doc, "<!DOCTYPE html>\n");
    fprintf(doc, "<html>\n");
    fprintf(doc, "<head>\n");
    fprintf(doc, "<meta charset='UTF-8'/>\n");
    fprintf(doc, "<meta name='viewport' content='width=device-width, initial-scale=1.0'/>\n");
    fprintf(doc, "<link href='../doc/mdb.min.css' rel='stylesheet'>\n");
    fprintf(doc, "<title>Schelling_MPI</title>\n");
    fprintf(doc, "<style>body {background-color: #1e2428}</style>\n");

    fprintf(doc, "</head>\n\n");

    fprintf(doc, "<body>");

    fprintf(doc, "<table id='dtVerticalScrollExample' class='table table-sm table-borderless' cellspacing='0' width='100%%'");

    for (i = 0; i < rows_size; i++) {
        fprintf(doc, "<tr>");
        for (j = 0; j < column_size; j++) {
            if (matrix[i * COLUMNS + j] == AGENT_O)
                fprintf(doc, "<td class='p-1' style='color:#ff4444'>%c</td>", matrix[i * COLUMNS + j]);
            else if (matrix[i * COLUMNS + j] == AGENT_X)
                fprintf(doc, "<td class='p-1' style='color:#0099CC'>%c</td>", matrix[i * COLUMNS + j]);
            else
                fprintf(doc, "<td class='p-1'>%c</td>", matrix[i * COLUMNS + j]);
        }
        fprintf(doc, "\n");
        fprintf(doc, "</tr>");
    }
    fprintf(doc, "</table>");
    fprintf(doc, "<script>\n");
    fprintf(doc, "$(document).ready(function () { $('#dtVerticalScrollExample').DataTable({ 'scrollY': '200px', 'scrollCollapse': true,}); $('.dataTables_length').addClass('bs-select'); });\n");
    fprintf(doc, "</script>\n");
    fprintf(doc, "</body>\n");
    fprintf(doc, "</html>");
    fclose(doc);
}

void err_finish(int *sendcounts, int *displacements, int *rows_per_process) {
    free(sendcounts);
    free(displacements);
    free(rows_per_process);
    MPI_Abort(MPI_COMM_WORLD, MPI_ERR_COUNT);
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