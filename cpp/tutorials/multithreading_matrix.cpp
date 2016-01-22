#include <thread>
#include <iostream>
#include <ctime>
#include <cstdlib>
#include <vector>

using namespace std;

void single_thread(int N);
void multi_thread(int N);
void code(int N, double &output, double &matrix, int start, int end);

int main()
{

  int X = 500;

  //time single_thread

  single_thread(X);
  multi_thread(X);


  return 0;
}


void single_thread(int N)
{
  double matrix [N][N];
  double output [N][N];

  //initialize
  for(int i=0; i < N; i++)
    for(int k=0; k < N; k++)
      output[i][k] = 0.0;

  //initialize random matrix
  srand( (unsigned)time( NULL ) );// set the seed
  for(int i=0; i < N; i++)
    for(int k=0; k < N; k++)
      matrix[i][k] = rand();

  clock_t begin = clock();
  for (int i = 0; i < N; i++)
  {
    for (int j = 0; j < N; j++)
    {
      for (int k = 0; k < N; k++)
      {
        output[i][j] += (matrix[i][k] * matrix[k][j]);
      }
    }
  }
  clock_t end = clock();
  cout << "Elapsed time for single threaded multiplication: " << double(end - begin) / CLOCKS_PER_SEC << endl;

  return;
}

void multi_thread(int N)
{
  int numthreads = thread::hardware_concurrency();

  int rows = N / numthreads;
  int extra = N % numthreads;
  int start = 0;
  int end = rows;
  double matrix [N][N];
  double output [N][N];

  //initialize random  & storage matrix
  srand( (unsigned)time( NULL ) );// set the seed
  for(int i=0; i < N; i++)
    for(int k=0; k < N; k++)
    {
      output[i][k] = 0.0;
      matrix[i][k] = rand();
    }

  clock_t begin = clock();

  auto code = [N,&output,&matrix](int start, int end) -> void
  {
    for (int i = start; i < end; i++)
      for (int j = 0; j < N; j++)
        for (int k = 0; k < N; k++)
          output[i][j] += (matrix[i][k] * matrix[k][j]);
  };


  vector<thread> workers;
  for (int t = 1; t <= numthreads; t++)
  {
    if (t == numthreads)
      end += extra;

    workers.push_back( thread(code, start, end) );
    start = end;
    end = start + rows;
  }

  for (thread& t : workers)
    t.join();

  clock_t c_end = clock();
  cout << "Elapsed time for single threaded multiplication: " << double(c_end - begin) / CLOCKS_PER_SEC << endl;
  return;
}

/*
void code(int N, double &output, double &matrix,int start, int end)
{
    for (int i = start; i < end; i++)
      for (int j = 0; j < N; j++)
        for (int k = 0; k < N; k++)
          output[i][j] += (matrix[i][k] * matrix[k][j]);
}
*/
