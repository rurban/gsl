/* define how a problem is split recursively */
#define GSL_LINALG_SPLIT(n) ((n >= 16) ? ((n + 8) / 16) * 8 : n / 2)

/* matrix size for crossover to Level 2 algorithms */
#define CROSSOVER              24
#define CROSSOVER_CHOLESKY     CROSSOVER
#define CROSSOVER_INVTRI       CROSSOVER
#define CROSSOVER_TRIMULT      CROSSOVER
