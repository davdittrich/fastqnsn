# Post-Mortem: `fastqnsn` v0.1.0 to v0.5.0 

## 1. The Journey at a Glance

The development of `fastqnsn` from its initial v0.1.0 commit to the highly optimized v0.5.0 release represents a masterclass in performance engineering, algorithmic refinement, and architectural pivoting. The primary goal was clear: implement the Rousseeuw-Croux $Q_n$ and $S_n$ estimators with bit-exact accuracy against the gold-standard `robustbase` package, but at a fraction of the computational cost.

The final result is a package that achieves **7.3x speedup for $S_n$** and **5.8x speedup for $Q_n$** at $n=100,000$, while delivering over **half a million windows per second** for streaming micro-samples ($n=10$).

This journey was not linear. It required aggressive experimentation, throwing away "fast" code when it failed exactness constraints, and deeply analyzing hardware bottlenecks.

---

## 2. What Went Well

### Data-Driven Scaling Strategy (The Hybrid Architecture)
The most successful architectural decision was acknowledging that **one size does not fit all**. Parallelizing an algorithm always introduces overhead. 
* We successfully implemented a **tri-modal execution model**:
  * **Micro-scale ($n \le 300$)**: Bruteforce $O(n^2)$ for $Q_n$ mapped perfectly to modern CPU caches.
  * **Mid-scale ($n < 10000$)**: Single-threaded $O(n \log n)$ bypassing thread-spawning overhead.
  * **Macro-scale ($n \ge 10000$)**: `RcppParallel` with Intel TBB, yielding 2.3x speedups on 16 threads.
* **Insight:** Scaling requires profiling the *boundaries* between approaches, not just the asymptotes.

### Algorithmic Upgrades Over Brute Force
In v0.5.0, replacing `std::nth_element` with the **Floyd-Rivest** selection algorithm yielded a ~30% reduction in comparisons. 
* **Insight:** When you hit the ceiling of compiler flags and memory layout, you must return to computer science fundamentals. Better algorithms beat better hardware.

### Memory Architecture (Arena Allocation)
Moving from multiple discrete heap allocations (`std::unique_ptr<T[]>`) to a single pre-allocated memory arena for the $Q_n$ worker arrays reduced OS-level allocation contention.
* **Insight:** In high-throughput C++, memory allocation is often the silent killer of performance. Stack allocation (for micro-samples) and arena allocation (for macro-samples) are essential.

### Strict Verification Scaffolding
We built highly aggressive benchmark and verification scripts (`verify_akinshin.R`, `benchmark_*.R`) early on. By forcing ourselves to match `robustbase` output bit-for-bit (e.g., discovering the need for `(float)` truncation in the $Q_n$ `whimed` loop), we prevented "fast but wrong" code from entering the main branch.

---

## 3. What Did Not Go Well (Failures & Pivots)

### The Zig Experiment (The Sunk Cost)
In v0.2.0, we attempted to write the high-performance kernels in **Zig**, aiming for seamless C-ABI integration and aggressive cross-compilation. 
* **The Failure:** While fast, the interop between R, Rcpp, and Zig proved fragile. Memory management across the FFI boundary was complex, and ensuring bit-exact behavior with R's `NA` semantics was a constant battle.
* **The Pivot:** In v0.3.0, we ripped out the Zig codebase entirely and rewrote the kernels in pure modern **C++17**, leveraging `RcppParallel` and `BH` (Boost). 
* **The Lesson:** "Boring" contiguous tech stacks (R + Rcpp + C++) often win over fragmented "blazing fast" polyglot stacks because the boundary overhead (mental and computational) destroys the theoretical gains.

### Amdahl's Law and the $S_n$ Synchronization Spike
Initially, we tried to over-parallelize the $S_n$ algorithm for medium $N$. We discovered a massive performance regression where the multi-threaded version was *slower* than the serial version due to thread synchronization overhead (the "thread-spawning spike").
* **The Fix:** We had to implement strict grain sizes and curb parallelism to 4 threads for medium workloads, pushing the parallel threshold up to $n=10000$.
* **The Lesson:** Default thread schedulers (like TBB) evaluate work dynamically, but the overhead of that evaluation can exceed the work itself for $O(n \log n)$ algorithms on small $N$.

### Over-Engineering the `whimed` Pivot
In v0.5.0, we attempted to optimize the weighted median (`whimed`) trial loop by replacing the simple midpoint pivot with a "median-of-3" approach to improve convergence.
* **The Reality:** Profiling showed the median-of-3 added branching overhead that actually *slowed down* execution, while offering zero tangible reduction in iterations. We reverted the change immediately.
* **The Lesson:** Always A/B test algorithmic "improvements." Theoretical convergence guarantees do not always translate to fewer CPU cycles, especially when branch prediction is involved.

---

## 4. Key Insights for the Future

1. **A/B Testing is Non-Negotiable in Systems Code:** 
   We successfully used a scientific approach to optimization: implement, benchmark, measure against baseline, and revert if the ROI is negative. We discarded branchless counting and software pipelining because the baseline `__builtin_expect` and sequential loops were heavily optimized by GCC already.
2. **Precision is an Algorithmic Constraint:** 
   The $Q_n$ algorithm is highly sensitive to floating-point truncation. Matching a reference implementation exactly (Akinshin 2022 / robustbase) required intentional precision "downgrades" (e.g., casting doubles to floats before comparison). You cannot optimize math without first establishing the exact precision semantics required.
3. **The Limits of Strong Scaling:**
   We hit an Amdahl's Law wall at around 2.3x speedup on 16 threads for $Q_n$. The serial sorting phase and the inherently sequential O(log n) `whimed` convergence loop dictate the maximum theoretical speedup. 
   * **Future work** must focus on completely parallelizing the sort (using GPU or better multi-core radix sorts) or finding a non-iterative approach to the trial evaluation.

## 5. Conclusion

The jump from v0.1.0 to v0.5.0 transformed `fastqnsn` from a functional prototype into a production-grade, state-of-the-art package. We achieved our goals not through a single "silver bullet," but through the compounded interest of dozens of micro-optimizations: hybrid thresholds, Boost Spreadsort, Floyd-Rivest selection, arena allocation, and strict compiler hints. 

The biggest meta-lesson? **The fastest code is the code that doesn't run.** (e.g., hoisting `std::max` out of loops, bypassing thread creation for small $N$, and using $O(n^2)$ brute-force to skip complex tree traversals entirely).
