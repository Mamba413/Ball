#include <algorithm>

bool sort_pair_compare(std::pair<double, int> a, std::pair<double, int> b)
{
  return a.first < b.first;
}
