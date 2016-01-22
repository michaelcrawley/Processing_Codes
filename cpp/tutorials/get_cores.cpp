#include <iostream>
#include <thread>

using namespace std;

int main()
{
  cout << "**Number of detected cores: " << thread::hardware_concurrency() << endl;
  return 0;
}
