#include <iostream>
#include <typeinfo>
#include <type_traits>
using namespace std;

#include "static_poly/static_poly_traits.hpp"
#include "static_poly/encap.hpp"
#include "static_poly/encap_manager.hpp"

template<
  class EncapManagerT1,
  class EncapManagerT2 = EncapManagerT1,
  class EncapManagerTraits1 = traits::encapmanager_traits<EncapManagerT1>,
  class EncapManagerTraits2 = traits::encapmanager_traits<EncapManagerT2>
>
class EncapManagerTransfer
{
  static_assert(is_base_of<EncapManager<typename EncapManagerTraits1::encap_type>, EncapManagerT1>::value,
                "first EncapManager of Transfer must be an EncapManager");
  static_assert(is_base_of<EncapManager<typename EncapManagerTraits2::encap_type>, EncapManagerT2>::value,
                "second EncapManager of Transfer must be an EncapManager");

  public:
    typedef EncapManagerTraits1 first_manager_traits;
    typedef EncapManagerTraits2 second_manager_traits;
    typedef typename first_manager_traits::encap_type first_encap_type;
    typedef typename second_manager_traits::encap_type second_encap_type;
    static_assert(is_convertible<first_encap_type, second_encap_type>::value,
                  "Encap of first Manager not convertible to Encap of second Manager");
    static_assert(is_convertible<second_encap_type, first_encap_type>::value,
                  "Encap of second Manager not convertible to Encap of first Manager");

    typedef typename first_manager_traits::time_type first_time_type;
    typedef typename first_manager_traits::spatial_type first_spatial_type;
    typedef typename second_manager_traits::time_type second_time_type;
    typedef typename second_manager_traits::spatial_type second_spatial_type;
};


int main(int argn, char** argv)
{
  Encapsulation<double, float> d;
  d.foo();
  cout << endl;

  ExtendedEncapsulation<double, int> e = d;
  e.foo();

  EncapManager<Encapsulation<float, char>> m;

  ExtendedEncapManager<float, double> e2;

  EncapManagerTransfer<EncapManager<Encapsulation<double, float>>> t1;
  EncapManagerTransfer<ExtendedEncapManager<double, float>> t2;
  EncapManagerTransfer<EncapManager<Encapsulation<double, float>>, ExtendedEncapManager<double, float>> t3;
}
