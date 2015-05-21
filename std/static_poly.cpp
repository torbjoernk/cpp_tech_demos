#include <iostream>
#include <typeinfo>
#include <type_traits>
using namespace std;


namespace traits
{
  template<class T>
  using is_temporal_type = is_floating_point<T>;
}


template<
  typename TimeT,
  typename SpacialT
>
class Encapsulation
{
  static_assert(traits::is_temporal_type<TimeT>::value,
                "temporal type must be a floating point");

  public:
    typedef TimeT    time_type;
    typedef SpacialT spacial_type;
    
    Encapsulation()
    {
      cout << "Encapsulation<" << typeid(time_type).name() << ", "
                               << typeid(spacial_type).name()
           << ">: empty ctor" << endl;
    }
    
    virtual spacial_type foo()
    {
      cout << "Encapsulation<" << typeid(time_type).name() << ", "
                               << typeid(spacial_type).name()
           << ">:foo()" << endl;
      return spacial_type(42);
    }
};

static_assert(
  is_same<
    double,
    Encapsulation<double, float>::time_type
  >::value, "");


template<
  typename TimeT,
  typename SpacialT
>
class ExtendedEncapsulation
  : public Encapsulation<TimeT, SpacialT>
{
  public:
    using typename Encapsulation<TimeT, SpacialT>::time_type;
    using typename Encapsulation<TimeT, SpacialT>::spacial_type;
    
    ExtendedEncapsulation()
    {
      cout << "ExtendedEncapsulation<" << typeid(time_type).name() << ", "
                                       << typeid(spacial_type).name()
           << ">: empty ctor" << endl;
    };
    
    template<typename OtherTimeT, typename OtherSpacialT>
    ExtendedEncapsulation(const Encapsulation<OtherTimeT, OtherSpacialT>&)
    {
      cout << "ExtendedEncapsulation<" << typeid(time_type).name() << ", "
                                       << typeid(spacial_type).name()
           << ">: copy ctor" << endl;
    };
    
    virtual spacial_type foo()
    {
      cout <<"ExtendedEncapsulation<" << typeid(time_type).name() << ", "
                                      << typeid(spacial_type).name()
           << ">:foo()" << endl;
      return spacial_type(21);
    }
};

static_assert(
  is_same<
    Encapsulation<double, float>::time_type,
    ExtendedEncapsulation<double, int>::time_type
  >::value, "");
static_assert(
  is_same<
    int,
    ExtendedEncapsulation<double, int>::spacial_type
  >::value, "");


template<
  template<
    typename...
  > class EncapT,
  typename... EncapTTs
>
class EncapManager
{
  static_assert(is_base_of<Encapsulation<EncapTTs...>, EncapT<EncapTTs...>>::value,
                "EncapManager requires an Encapsulation");

  public:
    typedef EncapT<EncapTTs...> encap_type;
};



int main(int argn, char** argv)
{
  Encapsulation<double, float> d;
  d.foo();
  cout << endl;
  
  ExtendedEncapsulation<double, int> e = d;
  e.foo();
  
  EncapManager<Encapsulation, size_t, char> m;
}
