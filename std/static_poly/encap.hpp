#ifndef _STATIC_POLY__ENCAP_HPP_
#define _STATIC_POLY__ENCAP_HPP_

#include <iostream>
#include <type_traits>
using namespace std;

#include "static_poly_traits.hpp"


template<
  typename TimeT,
  typename SpatialT
>
class Encapsulation
{
  static_assert(is_floating_point<TimeT>::value,
                "temporal type must be a floating point");
  static_assert(is_arithmetic<SpatialT>::value,
                "spatial type must be an arithmetic type");

  public:
    typedef TimeT    time_type;
    typedef SpatialT spatial_type;

    Encapsulation();
    virtual spatial_type foo();
};


template<
  typename TimeT,
  typename SpatialT
>
Encapsulation<TimeT, SpatialT>::Encapsulation()
{
  cout << "Encapsulation<" << typeid(time_type).name() << ", "
                           << typeid(spatial_type).name()
       << ">: empty ctor" << endl;
}

template<
  typename TimeT,
  typename SpatialT
>
typename Encapsulation<TimeT, SpatialT>::spatial_type
Encapsulation<TimeT, SpatialT>::foo()
{
  cout << "Encapsulation<" << typeid(time_type).name() << ", "
                           << typeid(spatial_type).name()
       << ">:foo()" << endl;
  return spatial_type(42);
}


static_assert(is_same<double, Encapsulation<double, float>::time_type>::value,
              "time type of Encapsulation is not set correctly");
static_assert(is_same<float, Encapsulation<double, float>::spatial_type>::value,
              "spatial type of Encapsulation is not set correctly");


template<
  typename TimeT,
  typename SpatialT
>
class ExtendedEncapsulation
  : public Encapsulation<TimeT, SpatialT>
{
  private:
    typedef Encapsulation<TimeT, SpatialT> _base_type;

  public:
    using typename   _base_type::time_type;
    using typename   _base_type::spatial_type;

    ExtendedEncapsulation();

    template<typename OtherTimeT, typename OtherSpatialT>
    ExtendedEncapsulation(const Encapsulation<OtherTimeT, OtherSpatialT>&);

    virtual spatial_type foo();
};

template<
  typename TimeT,
  typename SpatialT
>
ExtendedEncapsulation<TimeT, SpatialT>::ExtendedEncapsulation()
{
  cout << "ExtendedEncapsulation<" << typeid(time_type).name() << ", "
                                   << typeid(spatial_type).name()
       << ">: empty ctor" << endl;
};

template<
  typename TimeT,
  typename SpatialT
>
template<
  typename OtherTimeT, typename OtherSpatialT
>
ExtendedEncapsulation<TimeT, SpatialT>::ExtendedEncapsulation(const Encapsulation<OtherTimeT, OtherSpatialT>&)
{
  cout << "ExtendedEncapsulation<" << typeid(time_type).name() << ", "
                                   << typeid(spatial_type).name()
       << ">: copy ctor" << endl;
};

template<
  typename TimeT,
  typename SpatialT
>
typename ExtendedEncapsulation<TimeT, SpatialT>::spatial_type
ExtendedEncapsulation<TimeT, SpatialT>::foo()
{
  cout <<"ExtendedEncapsulation<" << typeid(time_type).name() << ", "
                                  << typeid(spatial_type).name()
       << ">:foo()" << endl;
  return spatial_type(21);
}

static_assert(is_same<Encapsulation<double, float>::time_type, ExtendedEncapsulation<double, int>::time_type>::value,
              "time type of ExtendedEncapsulation is not derived from Encapsulation");
static_assert(is_same<Encapsulation<double, int>::spatial_type, ExtendedEncapsulation<double, int>::spatial_type>::value,
              "spatial type of ExtendedEncapsulation is not derived from Encapsulation");


static_assert(traits::is_encapsulation<Encapsulation<double, int>>::value,
              "is_encapsulation trait is broken for base Encapsulation");
static_assert(traits::is_encapsulation<ExtendedEncapsulation<double, int>>::value,
              "is_encapsulation trait is broken for derived Encapsulations");


static_assert(is_same<traits::encapsulation_traits<Encapsulation<double, int>>::type, Encapsulation<double, int>>::value,
              "encapsulation_traits is broken for dependent type");
static_assert(is_same<traits::encapsulation_traits<Encapsulation<double, int>>::time_type, double>::value,
              "encapsulation_traits is broken for dependent time type");
static_assert(is_same<traits::encapsulation_traits<Encapsulation<double, int>>::spatial_type, int>::value,
              "encapsulation_traits is broken for dependent spatial type");

#endif // _STATIC_POLY__ENCAP_HPP_
