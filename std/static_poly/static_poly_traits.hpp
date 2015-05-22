#ifndef _STATIC_POLY__TRAITS_HPP_
#define _STATIC_POLY__TRAITS_HPP_

#include <type_traits>
using namespace std;


template<typename TimeT, typename SpatialT>
class Encapsulation;


namespace traits
{
  template<class EncapT>
  struct is_encapsulation
    : public integral_constant<bool, is_base_of<Encapsulation<typename EncapT::time_type, typename EncapT::spatial_type>, EncapT>::value>
  {};

  template<class EncapT, class Enabled=void>
  struct encapsulation_traits {};

  template<class EncapT>
  struct encapsulation_traits<EncapT, typename enable_if<is_encapsulation<EncapT>::value>::type> {
    typedef          EncapT               type;
    typedef typename EncapT::time_type    time_type;
    typedef typename EncapT::spatial_type spatial_type;
  };
}


template<class EncapT, class EncapTraitsT = traits::encapsulation_traits<EncapT>>
class EncapManager;


namespace traits
{
  template<class EncapManagerT, class Enabled=void>
  struct encapmanager_traits {};

  template<class EncapManagerT>
  struct is_encapmanager
    : public integral_constant<bool, is_base_of<EncapManager<typename EncapManagerT::encap_traits::type>, EncapManagerT>::value>
  {};

  template<class EncapManagerT>
  struct encapmanager_traits<EncapManagerT, typename enable_if<is_encapmanager<EncapManagerT>::value>::type>
  {
    typedef          EncapManagerT                            manager_type;
    typedef typename manager_type::encap_traits::type         encap_type;
    typedef typename manager_type::encap_traits::time_type    time_type;
    typedef typename manager_type::encap_traits::spatial_type spatial_type;
  };
}

#endif  // _STATIC_POLY__TRAITS_HPP_
