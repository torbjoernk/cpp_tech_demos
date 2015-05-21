#include <iostream>
using namespace std;


template<
  typename T1,
  typename... Ts
>
struct MixinOne
{
  T1 _m1;
};

template<
  typename T2,
  typename... Ts
>
struct MixinTwo
{
  T2 _m2;
};




int main(int argc, char** argv)
{
  
}
