// Private audit executable: defining Fermi integral; no production dependency.
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <boost/math/constants/constants.hpp>
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <array>
using MP=boost::multiprecision::cpp_dec_float_100;
std::array<MP,4> integral(const MP& eta) {
 if(eta < -150 || eta > 150) throw std::runtime_error("audit eta domain");
 const MP upper=sqrt(MP(400)+abs(eta));
 const MP norm=4/sqrt(boost::math::constants::pi<MP>());
 auto value=[&](const MP& t)->MP {const MP e=exp(eta-t*t);return MP(norm*t*t*e/(1+e));};
 auto deriv=[&](const MP& t)->MP {const MP e=exp(eta-t*t);return MP(norm*t*t*e/((1+e)*(1+e)));};
 using Q=boost::math::quadrature::gauss_kronrod<MP,61>;
 MP ev,ed;
 MP v=Q::integrate(value,MP(0),upper,15,MP("1e-40"),&ev);
 MP d=Q::integrate(deriv,MP(0),upper,15,MP("1e-40"),&ed);
 if(ev>abs(v)*MP("1e-32")||ed>abs(d)*MP("1e-32"))throw std::runtime_error("quadrature error estimate too large");
 return {v,d,ev,ed};
}
int main(){
 std::cout<<std::setprecision(60);char mode;MP x;
 while(std::cin>>mode>>x){
  try {
   if(mode=='f'){auto [v,d,ev,ed]=integral(x);std::cout<<x<<' '<<v<<' '<<d<<' '<<ev<<' '<<ed<<'\n';continue;}
   if(mode!='r'||x<=0||x<MP("1e-40")||x>100)throw std::runtime_error("audit ratio domain");
   MP eta=x>1?pow(3*sqrt(boost::math::constants::pi<MP>())*x/4,MP(2)/3):log(x);
   bool converged=false;
   for(int k=0;k<30;++k){auto [v,d,ev,ed]=integral(eta);MP delta=(v-x)/d;if(abs(delta)>10)delta=delta>0?MP(10):MP(-10);eta-=delta;if(abs(delta)<MP("1e-35")){converged=true;break;}}
   if(!converged)throw std::runtime_error("inverse did not converge");
   auto [v,d,ev,ed]=integral(eta);std::cout<<eta<<' '<<v<<' '<<d<<' '<<ev<<' '<<ed<<'\n';
  }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 2;}
 }
}
