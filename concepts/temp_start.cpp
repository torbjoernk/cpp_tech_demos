#include <iostream>
#include <typeinfo>
#include <type_traits>
using namespace std;


namespace mixins
{
  template<typename ValueT>
  struct NotCommunicateable
  {
    typedef ValueT comm_value_type;
    virtual void send() { cout << "NotCommunicatable::send()" << endl; }
    virtual void recv() { cout << "NotCommunicatable::recv()" << endl; }
  };

  template<typename ValueT>
  struct Communicateable
  {
    typedef ValueT comm_value_type;
    virtual void send() { cout << "Communicatable::send()" << endl; }
    virtual void recv() { cout << "Communicatable::recv()" << endl; }
  };

  template<typename ValueT>
  struct DerivedCommunicateable
    : Communicateable<ValueT>
  {
    virtual void send() { cout << "SpecialCommunicatable::send()" << endl; }
    virtual void recv() { cout << "SpecialCommunicatable::recv()" << endl; }
  };
}


template<
  typename ValueT,
  template<
    typename
  > class CommunicateableT = mixins::NotCommunicateable,
  typename... Ts
>
class Encapsulation
  : public CommunicateableT<ValueT>
{
  public:
    typedef ValueT value_type;
    typedef CommunicateableT<ValueT> comm_type;
};


template<
  typename ValueT,
  typename... Ts
>
class SpecialEncapsulation
  : public Encapsulation<ValueT, mixins::DerivedCommunicateable, Ts...>
{};



namespace traits
{
  template<typename>
  struct is_communicateable
  {
    static const bool value = false;
    typedef integral_constant<bool, value> type;
  };

  template<typename ValueT>
  struct is_communicateable<mixins::Communicateable<ValueT>>
  {
    static const bool value = true;
    typedef integral_constant<bool, value> type;
  };
  
  template<class CommT>
  struct communicateable_trait
  {
    using comm_type = typename CommT::comm_type;
    using valid_comm_trait = is_communicateable<CommT>;
  };

  template<class EncapT>
  struct encapsulation_trait
  {
    using value_type = typename EncapT::value_type;
    using comm_type = typename EncapT::comm_type;
  };
}


template<
  class EncapT
>
class EncapManager
{
  typedef typename traits::encapsulation_trait<EncapT> encap_traits;
};



int main(int argn, char** argv)
{
  Encapsulation<double> default_comm;
  default_comm.send();
  default_comm.recv();

  Encapsulation<double, mixins::Communicateable> custom_comm;
  custom_comm.send();
  custom_comm.recv();

  SpecialEncapsulation<double> special;
  special.send();
  special.recv();
  
  EncapManager<Encapsulation<double>> m;
}
