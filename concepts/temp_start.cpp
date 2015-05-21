#include <iostream>
#include <typeinfo>
#include <type_traits>
using namespace std;


namespace mixins
{
  template<typename>
  struct NotCommunicateable
  {
    typedef void value_type;
    virtual void send() { cout << "NotCommunicatable::send()" << endl; }
    virtual void recv() { cout << "NotCommunicatable::recv()" << endl; }
  };

  template<
    typename ValueT
  >
  struct Communicateable
  {
    typedef ValueT value_type;
    virtual void send() { cout << "Communicatable::send()" << endl; }
    virtual void recv() { cout << "Communicatable::recv()" << endl; }
  };

  template<
    typename ValueT
  >
  struct SpecialCommunicateable
    : Communicateable<ValueT>
  {
    typedef ValueT value_type;
    virtual void send() { cout << "SpecialCommunicatable::send()" << endl; }
    virtual void recv() { cout << "SpecialCommunicatable::recv()" << endl; }
  };
}


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
    using value_type = typename CommT::value_type;
    using valid_comm_trait = is_communicateable<CommT>;
  };
}


static_assert(is_same<double, typename traits::communicateable_trait<mixins::Communicateable<double>>::value_type>::value,
              "");


template<
  typename ValueT,
  template<
    typename
  > class CommunicateableT = mixins::NotCommunicateable,
  typename... Ts
>
class Encapsulation
  : public CommunicateableT<ValueT>
{};


template<
  typename ValueT,
  typename... Ts
>
class SpecialEncapsulation
  : public Encapsulation<ValueT, mixins::SpecialCommunicateable, Ts...>
{};


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
}
