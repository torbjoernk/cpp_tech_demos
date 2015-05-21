#include <iostream>
#include <string>
#include <typeinfo>
using namespace std;


template<
  typename T1,
  typename ...Ts
>
struct TemplateOne
{
  T1 x;
};

template<
  typename T1,
  typename T2,
  typename... Ts
>
struct TemplateTwo
  : public TemplateOne<T1, T2, Ts...>
{
  T2 y;

  TemplateTwo() = default;
  template<
    template<typename...> class OtherT,
    typename... OtherTs
  >
  TemplateTwo(const OtherT<OtherTs...>& other)
  {}
};


template<
  template<typename...> class DataT,
  typename T1,
  typename ...Ts
>
struct FactoryOne
{
  typedef DataT<T1, Ts...> data_type;
  data_type create() { return data_type(); }
};

template<
  typename T1,
  typename ...Ts
>
struct FactoryOne<TemplateTwo, T1, Ts...>
{
  typedef TemplateTwo<T1, Ts...> data_type;
  typedef true_type specialized;
  data_type create() { cout << "spezialized create" << endl; return data_type(); }
};



int main(int argc, char** argv)
{
  auto d1 = TemplateOne<double>();
  cout << typeid(d1).name() << " -- " << typeid(d1.x).name() << endl;

  auto d2 = TemplateTwo<double, int>();
  cout << typeid(d2).name() << " -- " << typeid(d2.x).name() << " -- " << typeid(d2.y).name() << endl;

  TemplateTwo<double, int> d1_cast = TemplateTwo<double, int>(d1);
  cout << typeid(d1_cast).name() << " -- " << typeid(d1_cast.x).name() << " -- " << typeid(d1_cast.y).name() << endl;

  auto fac1 = FactoryOne<TemplateOne, double>();
  cout << typeid(fac1).name() << " -- " << typeid(FactoryOne<TemplateOne, double>::data_type).name() << endl;

  auto fac1d = fac1.create();
  cout << typeid(fac1d).name() << " -- " << typeid(fac1d.x).name() << endl;

  auto fac2 = FactoryOne<TemplateTwo, double, int>();
  cout << typeid(fac2).name() << " -- " << typeid(FactoryOne<TemplateTwo, double, int>::data_type).name() << endl;
  static_assert(FactoryOne<TemplateTwo, double, int>::specialized::value, "did not use the template specialization");

  auto fac2d = fac2.create();
  cout << typeid(fac2d).name() << " -- " << typeid(fac2d.x).name() << " -- " << typeid(fac2d.y).name() << endl;
}
