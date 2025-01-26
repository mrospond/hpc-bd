#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "grower.h"
// #include "beehive.h"
#include <string.h>
#include <sys/time.h>
// #include "glider.h"

typedef struct {
    int **array;
    int numOfRows;
    int numOfColumns;
} board;

board rememberPreviousBoard[2];
board *currentBoardGlobal;
int generation = 0;

#define BIDIM 2
static int sizes[BIDIM];
static int subsizes[BIDIM];
static int starts[BIDIM];
static MPI_Datatype type, resizedtype;

#define SINGLE_PROCESS 1
#define ROOT 0
#define GHOST_CELLS 2
#define START_ROW 1500
#define START_COLUMN 1500

// void copyInitialPattern(board *boardGlobal){
//     int **board = boardGlobal->array;

//     for (int i = 1; i <= boardGlobal->numOfRows; i++) {
//         for (int j = 1; j <= boardGlobal->numOfColumns; j++) {
//             board[i][j] = 0;
//         }
//     }

//     for (int i = 0; i < GLIDER_HEIGHT; i++) {
//         for (int j = 0; j < GLIDER_WIDTH; j++) {
//             int targetRow = START_ROW + i;
//             int targetCol = START_COLUMN + j;

//             if (targetRow >= 0 && targetRow < boardGlobal->numOfRows &&
//                 targetCol >= 0 && targetCol < boardGlobal->numOfColumns) {
//                 board[targetRow][targetCol] = glider[i][j];
//             }
//         }
//     }
// }

void copyInitialPattern(board *boardGlobal){
    int **board = boardGlobal->array;

    for (int i = 1; i <= boardGlobal->numOfRows; i++) {
        for (int j = 1; j <= boardGlobal->numOfColumns; j++) {
            board[i][j] = 0;
        }
    }

    for (int i = 0; i < GROWER_HEIGHT; i++) {
        for (int j = 0; j < GROWER_WIDTH; j++) {
            int targetRow = START_ROW + i;
            int targetCol = START_COLUMN + j;

            if (targetRow >= 0 && targetRow < boardGlobal->numOfRows &&
                targetCol >= 0 && targetCol < boardGlobal->numOfColumns) {
                board[targetRow][targetCol] = grower[i][j];
            }
        }
    }
}

// void copyInitialPattern(board *boardGlobal){
//     int **board = boardGlobal->array;

//     for (int i = 1; i <= boardGlobal->numOfRows; i++) {
//         for (int j = 1; j <= boardGlobal->numOfColumns; j++) {
//             board[i][j] = 0;
//         }
//     }

//     for (int i = 0; i < BEEHIVE_HEIGHT; i++) {
//         for (int j = 0; j < BEEHIVE_WIDTH; j++) {
//             int targetRow = START_ROW + i;
//             int targetCol = START_COLUMN + j;

//             if (targetRow >= 0 && targetRow < boardGlobal->numOfRows &&
//                 targetCol >= 0 && targetCol < boardGlobal->numOfColumns) {
//                 board[targetRow][targetCol] = beehive[i][j];
//             }
//         }
//     }
// }

static void
world_print(board *world)
{
    int **cells = world->array;
    int row, col;

    for (row = 1; row <= world->numOfRows; row++)
    {
        for (col = 1; col <= world->numOfColumns; col++)
        {
            if (cells[row][col])
            {
                printf("O");
            }
            else
            {
                printf(" ");
            }
        }
        printf("\n");
    }
}

int countLiveCells(board *b) {
    int **board = b->array;

    int sumCells = 0;
    for (int i = 1; i <= b->numOfRows; i++) {
        for (int j = 1; j <= b->numOfColumns; j++) {
            sumCells = sumCells + board[i][j];
        }
    }

    return sumCells;
}

void nextGeneration(board *old, board *new) {
    int **oldBoard = old->array;

    for (int i = 1; i <= new->numOfRows; i++) {
        for (int j = 1; j <= new->numOfColumns; j++) {
            int newCell = 0;

            int upperCell = oldBoard[i - 1][j];
            int rightCell = oldBoard[i][j + 1];
            int leftCell = oldBoard[i][j - 1];
            int lowerCell = oldBoard[i + 1][j];
            int upperLeftCell = oldBoard[i - 1][j - 1];
            int upperRightCell = oldBoard[i - 1][j + 1];
            int lowerLeftCell = oldBoard[i + 1][j - 1];
            int lowerRightCell = oldBoard[i + 1][j + 1];

            int sum = upperCell + lowerCell + leftCell + rightCell + upperLeftCell + upperRightCell + lowerLeftCell + lowerRightCell;

            if (sum == 3) {
                newCell = 1;
            } else if (sum == 2) {
                newCell = oldBoard[i][j];
            } else {
                newCell = 0;
            }

            new->array[i][j] = newCell;
        }
    }
}

double timeSecs(void) {
    struct timeval tv;

    if (gettimeofday(&tv, 0) != 0) {
        fprintf(stderr, "could not do timing\n");
        exit(1);
    }

    return tv.tv_sec + (tv.tv_usec / 1000000.0);
}

int** malloc2D(int numRows, int numColumns) {
    int** board = malloc(numRows * sizeof(int *));
    if (board == NULL) {
        exit(1);
    }

    board[0] = malloc(numRows * numColumns * sizeof(int));
    if (board[0] == NULL) {
        exit(1);
    }

    for (int i = 1; i < numRows; i++) {
        board[i] = board[0] + i * numColumns;
    }
    return board;
}

void calculateCounts(int numOfRows, int numOfProcesses, int *counts) {
    int totalRows = numOfRows;

    int rowsPerProcess = numOfRows / numOfProcesses;
    for (int i = 0; i < numOfProcesses; i++) {
        counts[i] = rowsPerProcess;
    }

    if (numOfRows % numOfProcesses == 0) {
        return;
    }

    totalRows -= rowsPerProcess * numOfProcesses;
    int ind = 0;
    while (totalRows > 0) {
        counts[ind++]++;
        totalRows--;
    }
}

void prepareDataDistribution(int numOfProcesses, int *counts, int *displs, int boardWidth, int boardHeight) {
    sizes[0] = boardWidth + GHOST_CELLS; 
    sizes[1] = boardHeight + GHOST_CELLS;

    subsizes[0] = 1;
    subsizes[1] = boardHeight + GHOST_CELLS;

    starts[0] = 0; 
    starts[1] = 0;

    MPI_Type_create_subarray(BIDIM, sizes, subsizes, starts, MPI_ORDER_C, MPI_INT, &type);
    MPI_Type_create_resized(type, 0, (boardHeight + GHOST_CELLS) * sizeof(int), &resizedtype);
    MPI_Type_commit(&resizedtype);

    int currRow = 1;
    for (int i = 0; i < numOfProcesses; i++) {
        displs[i] = currRow; 
        currRow += counts[i];
    }
}

void initializeGlobalBoard(int rank, int boardWidth, int boardHeight) {
    currentBoardGlobal = (board *)malloc(sizeof(board));

    currentBoardGlobal->numOfRows = boardWidth;
    currentBoardGlobal->numOfColumns = boardHeight;
    currentBoardGlobal->array = malloc2D(boardWidth + GHOST_CELLS, boardHeight + GHOST_CELLS);

    if (rank == ROOT) {
        copyInitialPattern(currentBoardGlobal);
    }
}

void exchangeRows(int rank, int size, board *currentBoardLocal, int boardHeight, int numOfRowsPerMachine) {
    // we exchange rows
    int prevNeig = (rank + size - 1) % size;
    int nextNeig = (rank + 1) % size;

    if (rank % 2 == 0) {
        // first send then rcv
        MPI_Send(
            currentBoardLocal->array[numOfRowsPerMachine],
            boardHeight + GHOST_CELLS,
            MPI_INT,
            nextNeig,
            ROOT,
            MPI_COMM_WORLD);

        MPI_Recv(
            currentBoardLocal->array[numOfRowsPerMachine + 1],
            boardHeight + GHOST_CELLS,
            MPI_INT,
            nextNeig,
            ROOT,
            MPI_COMM_WORLD,
            MPI_STATUS_IGNORE);

        MPI_Send(
            currentBoardLocal->array[1],
            boardHeight + GHOST_CELLS,
            MPI_INT,
            prevNeig,
            ROOT,
            MPI_COMM_WORLD);

        MPI_Recv(
            currentBoardLocal->array[0],
            boardHeight + 2,
            MPI_INT,
            prevNeig,
            ROOT,
            MPI_COMM_WORLD,
            MPI_STATUS_IGNORE);
    }
    else {
        MPI_Recv(
            currentBoardLocal->array[0],
            boardHeight + GHOST_CELLS,
            MPI_INT,
            prevNeig,
            ROOT,
            MPI_COMM_WORLD,
            MPI_STATUS_IGNORE);

        MPI_Send(
            currentBoardLocal->array[1],
            boardHeight + GHOST_CELLS,
            MPI_INT,
            prevNeig,
            ROOT,
            MPI_COMM_WORLD);

        MPI_Recv(
            currentBoardLocal->array[numOfRowsPerMachine + 1],
            boardHeight + GHOST_CELLS,
            MPI_INT,
            nextNeig,
            ROOT,
            MPI_COMM_WORLD,
            MPI_STATUS_IGNORE);

        MPI_Send(
            currentBoardLocal->array[numOfRowsPerMachine],
            boardHeight + GHOST_CELLS,
            MPI_INT,
            nextNeig,
            ROOT,
            MPI_COMM_WORLD);
    }
}

int main(int argc, char *argv[]) {
    int size, rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int totalGenerations = 5000;
    int boardWidth = 3000;
    int boardHeight = 3000;

    int *counts = malloc(size * sizeof(int));
    int *displs = malloc(size * sizeof(int));

    calculateCounts(boardWidth, size, counts);
    prepareDataDistribution(size, counts, displs, boardWidth, boardHeight);

    int numOfRowsPerMachine;
    if (size == SINGLE_PROCESS) {
        numOfRowsPerMachine = boardWidth;
    }
    else {
        numOfRowsPerMachine = counts[rank];
    }

    for (int i = 0; i < 2; i++) {
        rememberPreviousBoard[i].numOfRows = numOfRowsPerMachine;
        rememberPreviousBoard[i].numOfColumns = boardHeight;
        rememberPreviousBoard[i].array = malloc2D(numOfRowsPerMachine + GHOST_CELLS, boardHeight + GHOST_CELLS);
    }
    board *currentBoardLocal = &rememberPreviousBoard[generation % 2];
    initializeGlobalBoard(rank, boardWidth, boardHeight);

    // if (rank == ROOT)
    // {
    //     printf("\ninitial board:\n\n");
    //     world_print(currentBoardGlobal);
    // }

    // Data distribution at the start
    if (size == SINGLE_PROCESS) {
        currentBoardLocal = currentBoardGlobal;
    }
    else {
        MPI_Scatterv(currentBoardGlobal->array[0], counts, displs,
                     resizedtype,
                     currentBoardLocal->array[1], (counts[rank]) * (boardHeight + GHOST_CELLS), resizedtype,
                     0, MPI_COMM_WORLD);
    }

    // Collective operation may not have the effect of synchronising processes so need a barrier
    // to ensure that we measure time correctly
    MPI_Barrier(MPI_COMM_WORLD);
    double start = timeSecs();
    for (generation = 1; generation < totalGenerations; generation++) {
        board *nextBoardLocal;
        nextBoardLocal = &rememberPreviousBoard[generation % 2];

        if (size != SINGLE_PROCESS) {
            exchangeRows(rank, size, currentBoardLocal, boardHeight, numOfRowsPerMachine);
        }
        nextGeneration(currentBoardLocal, nextBoardLocal);

        // if (size == SINGLE_PROCESS)
        //     {
        //         currentBoardGlobal = currentBoardLocal;
        //     }
        //     else
        //     {
        //         MPI_Gatherv(currentBoardLocal->array[1],
        //                     (counts[rank]) * (boardWidth + GHOST_CELLS),
        //                     MPI_INT,
        //                     currentBoardGlobal->array[0],
        //                     counts,
        //                     displs,
        //                     resizedtype,
        //                     ROOT,
        //                     MPI_COMM_WORLD);
        //     }

        //     if (rank == ROOT)
        //     {
        //         printf("\nat time step %d:\n\n", generation);
        //         world_print(currentBoardGlobal);
        //     }

        currentBoardLocal = nextBoardLocal;
    }
    // We need a barrier here to measure the time it takes all processes to reach this point
    MPI_Barrier(MPI_COMM_WORLD);
    double end = timeSecs();

    // gather to sum live cells
    if (size == SINGLE_PROCESS) {
        currentBoardGlobal = currentBoardLocal;
    }
    else {
        MPI_Gatherv(currentBoardLocal->array[1],
                    (counts[rank]) * (boardHeight + GHOST_CELLS),
                    MPI_INT,
                    currentBoardGlobal->array[0],
                    counts,
                    displs,
                    resizedtype,
                    ROOT,
                    MPI_COMM_WORLD);
    }

    if (rank == ROOT) {
        printf("Population is %d cells\n", countLiveCells(currentBoardGlobal));
        fprintf(stderr, "%d workers and %d generations took %10.3f s\n", size, totalGenerations, end - start);
    }

    for (int i = 0; i < 2; i++) {
        free(rememberPreviousBoard[i].array[0]);
        free(rememberPreviousBoard[i].array);
    }

    MPI_Finalize();

    return 0;
}
