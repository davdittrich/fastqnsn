# fastqnsn 0.5.0

## Performance Improvements

* **Floyd-Rivest selection algorithm** replaces `std::nth_element` throughout, achieving ~30% fewer comparisons for Sn low-median and Qn final selection.
* **Arena memory allocation** for both Qn and Sn: reduces heap allocation overhead by consolidating multiple allocations into a single arena.
* **Cached `std::max` in Sn inner loop**: eliminates redundant double-evaluation of the min-diameter computation, yielding ~10% Sn improvement.
* **Raised naive Qn threshold to n=300**: the O(n²) brute-force path now covers more of the streaming micro-window regime, delivering 114k windows/sec at n=200.
* **`(float)` truncation in whimed**: matches robustbase precision semantics for bit-exact accuracy.
* **`__builtin_expect` branch hints** on hot counting loops in QnCountWorker and QnRefineWorker.

## Benchmark Results (vs robustbase, n=100,000)

| Estimator | robustbase | fastqnsn | Speedup |
|-----------|-----------|----------|---------|
| Sn        | 10.5 ms   | 1.45 ms  | **7.3x**|
| Qn        | 93.1 ms   | 16.0 ms  | **5.8x**|

## Throughput (windows per second)

| n   | Sn WPS  | Qn WPS  |
|-----|---------|---------|
| 10  | 683k    | 553k    |
| 50  | 541k    | 425k    |
| 200 | 283k    | 114k    |

# fastqnsn 0.4.0

* Initial optimized release with hybrid serial/parallel architecture.
* RcppParallel-based parallel counting and refinement for large N.
* Boost Spreadsort for medium N, TBB parallel sort for large N.
* Corrected consistency constants and Akinshin (2022) bias corrections.
