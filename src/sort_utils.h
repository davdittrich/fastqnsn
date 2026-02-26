#ifndef SORT_UTILS_H
#define SORT_UTILS_H

#include "thresholds.h"
#include <algorithm>
#include <boost/sort/spreadsort/spreadsort.hpp>
#include <tbb/parallel_sort.h>
#include <type_traits>

namespace fastqnsn {

template <typename Iterator> void optimized_sort(Iterator begin, Iterator end) {
  using T = typename std::iterator_traits<Iterator>::value_type;
  size_t n = std::distance(begin, end);

  if (n < 2)
    return;

  if constexpr (std::is_floating_point_v<T>) {
    if (n <= SORT_BOOST_THRESHOLD) {
      std::sort(begin, end);
    } else if (n < SORT_TBB_FLOAT_THRESHOLD) {
      boost::sort::spreadsort::float_sort(begin, end);
    } else {
      tbb::parallel_sort(begin, end);
    }
  } else if constexpr (std::is_integral_v<T>) {
    if (n <= SORT_BOOST_THRESHOLD) {
      std::sort(begin, end);
    } else if (n < SORT_TBB_INT_THRESHOLD) {
      boost::sort::spreadsort::integer_sort(begin, end);
    } else {
      tbb::parallel_sort(begin, end);
    }
  } else {
    if (n < SORT_TBB_FLOAT_THRESHOLD) {
      std::sort(begin, end);
    } else {
      tbb::parallel_sort(begin, end);
    }
  }
}

} // namespace fastqnsn

#endif
