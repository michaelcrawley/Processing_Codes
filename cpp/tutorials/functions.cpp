#include <iostream>
using namespace std;

// function declaration
void increment(int &num);

int main ()
{
   // local variable declaration:
   int x = 1;

   cout << "Original value is : " << x << endl;
   increment(x);
   cout << "Incremented value is : " << x << endl;

   return 0;
}

// function returning the max between two numbers
void increment(int &num)
{
   // local variable declaration
   num = num +1;


   return;
}
