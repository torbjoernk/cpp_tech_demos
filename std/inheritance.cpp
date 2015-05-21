#include <iostream>
#include <string>
using namespace std;


class AbstractBaseA {
  protected:
    string _member;

  public:
    AbstractBaseA()
    {
      cout << "AbstractBaseA::AbstractBaseA()" << endl
           << "  _member: '" << this->_member << "'" << endl;
    }

    explicit AbstractBaseA(const string& mem)
      : _member(mem)
    {
      cout << "AbstractBaseA::AbstractBaseA('" << mem << "')" << endl
           << "  _member: '" << this->_member << "'" << endl;
    }

    virtual ~AbstractBaseA()
    {
      cout << "AbstractBaseA::~AbstractBaseA()" << endl
           << "  _member: '" << this->_member << "'" << endl;
    }

    virtual void func() = 0;
};


class CommonBaseA
  : public AbstractBaseA
{
  public:
    CommonBaseA()
      : AbstractBaseA()
    {
      cout << "CommonBaseA::CommonBaseA()" << endl
           << "  _member: '" << this->_member << "'" << endl;
    }
    explicit CommonBaseA(const string& mem)
      : AbstractBaseA(mem)
    {
      cout << "CommonBaseA::CommonBaseA('" << mem << "')" << endl
           << "  _member: '" << this->_member << "'" << endl;
    }

    ~CommonBaseA()
    {
      cout << "CommonBaseA::~CommonBaseA()" << endl
           << "  _member: '" << this->_member << "'" << endl;
    }

    virtual void func()
    {
      cout << "CommonBaseA::func()" << endl
           << "  _member: '" << this->_member << "'" << endl;
    }
};


class DerivedA
  : public AbstractBaseA
{
  public:
    DerivedA()
      : AbstractBaseA()
    {
      cout << "DerivedA::DerivedA()" << endl
           << "  _member: '" << this->_member << "'" << endl;
    }
    explicit DerivedA(const string& mem)
      : AbstractBaseA(mem)
    {
      cout << "DerivedA::DerivedA('" << mem << "')" << endl
           << "  _member: '" << this->_member << "'" << endl;
    }

    virtual ~DerivedA()
    {
      cout << "DerivedA::~DerivedA()" << endl
           << "  _member: '" << this->_member << "'" << endl;
    }

    virtual void func()
    {
      cout << "DerivedA::func()" << endl
           << "  _member: '" << this->_member << "'" << endl;
    }
};

class DerivedB
  : public CommonBaseA
{
  public:
    DerivedB()
      : CommonBaseA()
    {
      cout << "DerivedB::DerivedB()" << endl
           << "  _member: '" << this->_member << "'" << endl;
    }
    explicit DerivedB(const string& mem)
      : CommonBaseA(mem)
    {
      cout << "DerivedB::DerivedB('" << mem << "')" << endl
           << "  _member: '" << this->_member << "'" << endl;
    }

    virtual ~DerivedB()
    {
      cout << "DerivedB::~DerivedB()" << endl
           << "  _member: '" << this->_member << "'" << endl;
    }
};


int main(const int argn, const char** argv)
{
  {
    cout << "---" << endl;
    CommonBaseA base_a;
    base_a.func();
  }
  {
    cout << "---" << endl;
    CommonBaseA base_a("test string");
    base_a.func();
  }

  {
    cout << "---" << endl;
    DerivedA derived_a;
    derived_a.func();
  }
  {
    cout << "---" << endl;
    DerivedA derived_a("test string");
    derived_a.func();
  }

  {
    cout << "---" << endl;
    DerivedB derived_b;
    derived_b.func();
  }
  {
    cout << "---" << endl;
    DerivedB derived_b("test string");
    derived_b.func();
  }
}
