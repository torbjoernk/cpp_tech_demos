#include <cstdlib>
#include <exception>
#include <stdexcept>
#include <string>
using namespace std;


class CustomError
  : public runtime_error
{
  protected:
    string msg;

  public:
    explicit CustomError(const string& msg)
      : runtime_error(msg)
    {
      this->msg = string("prefix: ") + string(runtime_error::what());
    }

    virtual const char* what() const throw()
    {
      return this->msg.c_str();
    }
};


int main(int argc, char** argv)
{
  throw CustomError("custom message");
  exit(0);
}
