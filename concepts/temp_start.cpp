#include <iostream>
#include <typeinfo>
#include <type_traits>
using namespace std;


namespace policies
{
  template<typename>
  class NotCommunicateable
  {
    public:
      typedef void value_type;

      virtual void send() {};
      virtual void recv() {};
  };

  template<
    typename ValueT
  >
  class Communicateable
  {
    public:
      typedef ValueT value_type;

      virtual void send() { cout << "Communicatable::send()" << endl; }
      virtual void recv() { cout << "Communicatable::recv()" << endl; }
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
  struct is_communicateable<policies::Communicateable<ValueT>>
  {
    static const bool value = true;
    typedef integral_constant<bool, value> type;
  };

  template<class CommT>
  struct communicateable_trait
  {
    using value_type = typename CommT::value_type;
    using valid_comm = typename is_communicateable<CommT>::type;
  };
}

static_assert(is_same<double, typename traits::communicateable_trait<policies::Communicateable<double>>::value_type>::value,
              "");


template<
  typename ValueT,
  template<
    typename
  > class CommunicatorT = policies::NotCommunicateable,
  typename... Ts
>
class Encapsulation
  : public CommunicatorT<ValueT>
{};


int main(int argn, char** argv)
{
  Encapsulation<double> default_comm;
  default_comm.send();
  default_comm.recv();

  Encapsulation<double, policies::Communicateable> custom_comm;
  custom_comm.send();
  custom_comm.recv();
}
