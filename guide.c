      |            |
------------------------------
      |            |
      |   P(i,j)   |
      |            |
------------------------------
      |            |

// am de trimis si primit de la caolturi 2x2 elemente

int CalculeazaRangulProcesului(int rang_i, int rang_j)
{
return MPI_....
}

void NextGeneration(int nr, int r, int nc, int c, MatriceInitiala, &MatriceFinala)
{
# define max(a, b) ((a > b) ? (a) : (b))
int DIMENSIUNE_MAXIMA = mnax (nr / r + 4, nc / c + 4);
char *Buffer1 = malloc(DIMENSIUNE_MAXIMA);
char *Buffer2 = malloc(DIMENSIUNE_MAXIMA);
rang_i = MPI....
rang_j = MPI.....
[ etapa 1. trimit borderul de latime 2 la vecini ]
// Daca rang_i != 0 (nu e stanga)
if (rang_i != 0)
  {
  // Buffer1 = Coloana 0
  for (int i = 0; i < nr / r + 4 - 1; i++)
    Buffer1[i] = MatriceInitiala[i][0];
  // Buffer2 = Coloana 1
  for (int i = 0; i < nr / r + 4 - 1; i++)
    Buffer2[i] = MatriceInitiala[i][1];
  // trimite doua coloane din stanga la P(i-1, j)
  int dest = CalculeazaRangulProcesului(i - 1, j)
  int MPI_Send(Buffer1, nr / r + 4, MPI_Char, dest, 0, MPI_Comm)
  int MPI_Send(Buffer2, nr / r + 4, MPI_Char, dest, 0, MPI_Comm)
  }
Daca j != 0 (nu e jos)
  Buffer1 = Lm-1
  Buffer2 = Lm-2
  // trimite doua randuri din jos la P(i, j-1)
  int dest = CalculeazaRangulProcesului(i, j - 1)
  int MPI_Send(Buffer1, nc / c + 4, MPI_Char, dest, 0, MPI_Comm)
  int MPI_Send(Buffer2, nc / c + 4, MPI_Char, dest, 0, MPI_Comm)
Daca rang_i != r - 1 (nu e dreapta)
  Buffer1 = Cn-1
  Buffer2 = Cn-2
  // trimite doua coloane din dreapta la P(i+1, j)
  int dest = CalculeazaRangulProcesului(i + 1, j)
  int MPI_Send(Buffer1, nr / r + 4, MPI_Char, dest, 0, MPI_Comm)
  int MPI_Send(Buffer2, nr / r + 4, MPI_Char, dest, 0, MPI_Comm)
Daca rang_j != c - 1 (nu e sus)
  Buffer1 = C0
  Buffer2 = C1
  // trimite doua randuri din sus la P(i, j+1)
  int dest = CalculeazaRangulProcesului(i, j + 1)
  int MPI_Send(Buffer1, nr / r + 4, MPI_Char, dest, 0, MPI_Comm)
  int MPI_Send(Buffer2, nr / r + 4, MPI_Char, dest, 0, MPI_Comm)
daca rang_i != 0 si rang_j != 0
  {
  // trimite la vecinul stanga_jos
  int dest = CalculeazaRangulProcesului(i - 1, j - 1);
  Buffer1[0] = MatriceInitiala[0][0];
  Buffer1[1] = MatriceInitiala[0][1]; 
  Buffer1[2] = MatriceInitiala[1][0];
  Buffer1[3] = MatriceInitiala[1][1];
  int MPI_Send(Buffer1, 4, MPI_Char, dest, 0, MPI_Comm)
  }
daca rang_i != 0 si rang_j != nc / c - 1
  trimite la vecinul stanga_jos
daca rang_i != nr / r - 1 si rang_j != 0
  trimite la vecinul stanga_jos
daca rang_i != nr / r - 1 si rang_j != nc / c - 1
  trimite la vecinul stanga_jos
[ etapa 2. primesc  borderul de latime 2 la vecini ]
Daca rang_i != 0 (nu e stanga)
  {
  // asteapta doua coloane din stanga de la P(i-1, j)
  int source = CalculeazaRangulProcesului(i - 1, j)
  int status;
  int MPI_Recv(Buffer1, nr / r + 4, MPI_Char, source, 0, MPI_Comm, status)
  int MPI_Recv(Buffer2, nr / r + 4, MPI_Char, source, 0, MPI_Comm, status)
  // Coloana 0 = Buffer1
  for (int i = 0; i < nr / r + 4 - 1; i++)
    MatriceInitiala[i][0] = Buffer1[i];
  // Coloana 1 = Buffer2
  for (int i = 0; i < nr / r + 4 - 1; i++)
    MatriceInitiala[i][1] = Buffer2[i];
  }
Daca rang_j != 0 (nu e jos)
  asteapta doua randuri din jos de la P(i, j-1)
  L1
  L2
Daca rang_i != r - 1 (nu e dreapta)
  asteapta doua coloane din dreapta de la P(i+1, j)
  C3
  C4
Daca rang_j != c - 1 (nu e sus)
  asteapta doua randuri din sus de la P(i, j+1)
  L3
  L4
daca rang_i != 0 si rang_j != 0
  {
  // primeste de la vecinul stanga_jos
  int sursa = CalculeazaRangulProcesului(i - 1, j - 1);
  int MPI_Recv(Buffer1, 4, MPI_Char, sursa, 0, MPI_Comm)
  MatriceInitiala[0][0] = Buffer1[0];
  MatriceInitiala[0][1] = Buffer1[1];
  MatriceInitiala[1][0] = Buffer1[2];
  MatriceInitiala[1][1] = Buffer1[3];
  }
daca rang_i != 0 si rang_j != nc / c - 1
  trimite la vecinul stanga_jos
daca rang_i != nr / r - 1 si rang_j != 0
  trimite la vecinul stanga_jos
daca rang_i != nr / r - 1 si rang_j != nc / c - 1
  trimite la vecinul stanga_jos
[ etapa 3. calculele mele mai putin pe boder]
for (int i = 2; i <= nr / r + 2; i++)   // 0, 1, nr/r + 4 - 1, nr/r + 4 - 2 sunt borduri
  for (int j = 2; j <= nc / c + 2; j++) // 0, 1, nc/c + 4 - 1, nc/c + 4 - 2 sunt borduri
    MatriceNou[i][j] = Stensil(nr / r - 1, nc / c - 1, MatriceaVeche, i, j);
free(Buffer1);
free(Buffer2);
}