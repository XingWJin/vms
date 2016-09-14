#include "gmodel.hpp"
#include <sstream>
#include <cassert>

int main(int argc, char** argv)
{
  assert(argc == 2);
  int elems = atoi(argv[1]);
  double size = 1.0/(double)elems;
  gmod::default_size = size;
  auto s = gmod::new_square(
      gmod::Vector{0,0,0},
      gmod::Vector{1,0,0},
      gmod::Vector{0,1,0});
  std::ostringstream geo_name;
  std::ostringstream dmg_name;
  geo_name << "square_" << elems << ".geo";
  dmg_name << "square_" << elems << ".dmg";
  gmod::write_closure_to_geo(s, geo_name.str().c_str());
  gmod::write_closure_to_dmg(s, dmg_name.str().c_str());
}
