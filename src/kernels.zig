const std = @import("std");

/// Sn Naive (N^2)
/// For each i, calculates the median of |x[i] - x[j]|
/// Results are written to out_ptr
export fn zig_sn_naive(x_ptr: [*]const f64, n: usize, out_ptr: [*]f64) void {
    const out = out_ptr[0..n];
    const target_rank = (n - 1) / 2;

    for (0..n) |i| {
        out[i] = zig_sn_select_k(x_ptr, n, i, target_rank + 1);
    }
}

/// Sn Selector (Kernel)
/// For a given i, calculates the median of |x[i] - x[j]| for j in [0, n-1]
/// x must be sorted ascending.
export fn zig_sn_select_k(x_ptr: [*]const f64, n: usize, i: usize, k: usize) f64 {
    const x = x_ptr[0..n];
    const val_i = x[i];

    // Branchless binary search for the k-th smallest |x[i] - x[j]|
    // Utilizing the property that |x[i] - x[j]| is unimodal (v-shaped) in j
    var low: f64 = 0.0;
    var high: f64 = @max(val_i - x[0], x[n - 1] - val_i);

    for (0..64) |_| {
        const mid = low + (high - low) / 2.0;
        // Count pairs |x[i] - x[j]| <= mid
        // This is equivalent to finding j range [L, R] such that x[j] is in [val_i - mid, val_i + mid]
        var count: usize = 0;

        // Binary search for bounds (or two-pointer if we were doing this for all i)
        // Since this is for a single i, we use std.sort.binarySearch
        const left_bound = val_i - mid;
        const right_bound = val_i + mid;

        var l: usize = 0;
        var r: usize = n;
        while (l < r) {
            const m = l + (r - l) / 2;
            if (x[m] < left_bound) {
                l = m + 1;
            } else {
                r = m;
            }
        }
        const start = l;

        l = 0;
        r = n;
        while (l < r) {
            const m = l + (r - l) / 2;
            if (x[m] <= right_bound) {
                l = m + 1;
            } else {
                r = m;
            }
        }
        const end = l;
        count = end - start;

        if (count >= k + 1) { // +1 because |x[i]-x[i]| = 0 is always included
            high = mid;
        } else {
            low = mid;
        }
    }
    return high;
}

/// Qn Counter (Upper triangle only)
/// Counts pairs (i, j) with 0 <= j < i < n such that x[i] - x[j] <= trial.
/// x must be sorted ascending.
export fn zig_qn_count(x_ptr: [*]const f64, n: usize, trial: f64) u64 {
    const x = x_ptr[0..n];
    var count: u64 = 0;
    var j: usize = 0;
    for (1..n) |i| {
        while (j < i and (x[i] - x[j]) > trial) {
            j += 1;
        }
        count += (i - j);
    }
    return count;
}

/// Weighted Median Helper for JM
fn whimed(a: []f64, iw: []i32, n: usize, target: i64) f64 {
    if (n == 0) return 0.0;
    if (n == 1) return a[0];

    var l: usize = 0;
    var r: usize = n - 1;
    var t = target;

    while (l < r) {
        const pivot = a[l + (r - l) / 2];
        var i = l;
        var j = l;
        while (j <= r) : (j += 1) {
            if (a[j] < pivot) {
                const tmp_a = a[i];
                a[i] = a[j];
                a[j] = tmp_a;
                const tmp_w = iw[i];
                iw[i] = iw[j];
                iw[j] = tmp_w;
                i += 1;
            }
        }
        var wleft: i64 = 0;
        for (l..i) |idx| wleft += iw[idx];

        if (wleft > t) {
            r = i - 1;
        } else {
            t -= wleft;
            // Check for elements equal to pivot
            var i_eq = i;
            var j_eq = i;
            while (j_eq <= r) : (j_eq += 1) {
                if (a[j_eq] == pivot) {
                    const tmp_a = a[i_eq];
                    a[i_eq] = a[j_eq];
                    a[j_eq] = tmp_a;
                    const tmp_w = iw[i_eq];
                    iw[i_eq] = iw[j_eq];
                    iw[j_eq] = tmp_w;
                    i_eq += 1;
                }
            }
            var weq: i64 = 0;
            for (i..i_eq) |idx| weq += iw[idx];
            if (wleft + weq > t) {
                return pivot;
            } else {
                t -= weq;
                l = i_eq;
            }
        }
    }
    return a[l];
}

/// Advanced Qn Selector (Johnson-Mizoguchi)
/// Returns the k-th smallest difference x[i] - x[j] for j < i.
export fn zig_qn_jm_select(x_ptr: [*]const f64, n: usize, k: u64, work_ptr: [*]f64, iweight_ptr: [*]i32, left_ptr: [*]i32, right_ptr: [*]i32) f64 {
    const x = x_ptr[0..n];
    const left = left_ptr[0..n];
    const right = right_ptr[0..n];
    const work = work_ptr[0..n];
    const iweight = iweight_ptr[0..n];

    for (0..n) |i| {
        left[i] = 1;
        right[i] = @intCast(i);
    }

    var nL: u64 = 0;
    var nR: u64 = @intCast(n * (n - 1) / 2);

    while (nR - nL > n) {
        var m: usize = 0;
        for (1..n) |i| {
            if (left[i] <= right[i]) {
                const w = right[i] - left[i] + 1;
                const jj = left[i] + @divFloor(w, 2);
                work[m] = x[i] - x[i - @as(usize, @intCast(jj))];
                iweight[m] = w;
                m += 1;
            }
        }

        const trial = whimed(work[0..m], iweight[0..m], m, @intCast(@divFloor(nR - nL, 2)));

        // Count sumP (elements < trial) and sumQ (elements <= trial)
        var sumP: u64 = 0;
        var sumQ: u64 = 0;
        
        // Use two-pointer for efficient O(n) counting
        var jp: usize = 0;
        var jq: usize = 0;
        for (1..n) |i| {
            while (jp < i and (x[i] - x[jp]) >= trial) jp += 1;
            sumP += (i - jp);
            while (jq < i and (x[i] - x[jq]) > trial) jq += 1;
            sumQ += (i - jq);
        }

        if (k <= sumP) {
            for (1..n) |i| {
                // Update right[i] to be the largest jj such that x[i] - x[i-jj] < trial
                // This means x[i-jj] > x[i] - trial
                var cur_jj = right[i];
                while (cur_jj >= left[i] and (x[i] - x[i - @as(usize, @intCast(cur_jj))]) >= trial) {
                    cur_jj -= 1;
                }
                right[i] = cur_jj;
            }
            nR = sumP;
        } else if (k > sumQ) {
            for (1..n) |i| {
                // Update left[i] to be the smallest jj such that x[i] - x[i-jj] > trial
                var cur_jj = left[i];
                while (cur_jj <= right[i] and (x[i] - x[i - @as(usize, @intCast(cur_jj))]) <= trial) {
                    cur_jj += 1;
                }
                left[i] = cur_jj;
            }
            nL = sumQ;
        } else {
            return trial;
        }
    }

    // Final selection on candidates
    var m_final: usize = 0;
    for (1..n) |i| {
        var jj = left[i];
        while (jj <= right[i]) : (jj += 1) {
            work[m_final] = x[i] - x[i - @as(usize, @intCast(jj))];
            m_final += 1;
        }
    }

    const final_slice = work[0..m_final];
    std.sort.pdq(f64, final_slice, {}, std.sort.asc(f64));
    return final_slice[k - nL - 1];
}
