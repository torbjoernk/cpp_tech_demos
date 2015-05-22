#ifndef _STATIC_POLY__ENCAP_MANAGER_HPP_
#define _STATIC_POLY__ENCAP_MANAGER_HPP_

#include <iostream>
#include <type_traits>
using namespace std;

#include "static_poly_traits.hpp"
#include "encap.hpp"


template<
  class EncapT,
  class EncapTraitsT
>
class EncapManager
{
  public:
    typedef EncapTraitsT encap_traits;
};

static_assert(is_same<Encapsulation<double, float>::time_type, EncapManager<Encapsulation<double, float>>::encap_traits::time_type>::value,
              "time type of EncapManager is not extracted from Encapsulation's time type");


template<
  typename TimeT,
  typename SpatialT
>
class ExtendedEncapManager
  : public EncapManager<ExtendedEncapsulation<TimeT, SpatialT>>
{
  typedef EncapManager<ExtendedEncapsulation<TimeT, SpatialT>> _base_type;

  public:
    using typename _base_type::encap_traits;
};



#endif  // _STATIC_POLY__ENCAP_MANAGER_HPP_
