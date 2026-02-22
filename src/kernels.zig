const std = @import("std");

/// Sn Naive (N^2)
/// For each i, calculates the median of |x[i] - x[j]|
/// Results are written to out_ptr
export fn zig_sn_naive(x_ptr: [*]const f64, n: usize, out_ptr: [*]f64) void {
    const out = out_ptr[0..n];
    const target_rank = n / 2 + 1;

    for (0..n) |i| {
        out[i] = zig_sn_select_k(x_ptr, n, i, target_rank);
    }
}

/// Sn Selector (Kernel)
/// For a given i, calculates the median of |x[i] - x[j]| for j in [0, n-1]
/// x must be sorted ascending.
export fn zig_sn_select_k(x_ptr: [*]const f64, n: usize, i: usize, k: usize) f64 {
    const x = x_ptr[0..n];
    const val_i = x[i];

    var low: f64 = 0.0;
    var high: f64 = @max(val_i - x[0], x[n - 1] - val_i);

    for (0..64) |_| {
        const mid = low + (high - low) / 2.0;
        const left_bound = val_i - mid;
        const right_bound = val_i + mid;

        var l_idx: usize = 0;
        var r_idx: usize = n;
        while (l_idx < r_idx) {
            const m = l_idx + (r_idx - l_idx) / 2;
            if (x[m] < left_bound) {
                l_idx = m + 1;
            } else {
                r_idx = m;
            }
        }
        const start = l_idx;

        l_idx = 0;
        r_idx = n;
        while (l_idx < r_idx) {
            const m = l_idx + (r_idx - l_idx) / 2;
            if (x[m] <= right_bound) {
                l_idx = m + 1;
            } else {
                r_idx = m;
            }
        }
        const end = l_idx;
        const count = end - start;

        if (count >= k) {
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
            r = @max(l, i - 1);
        } else {
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
                t -= (wleft + weq);
                l = i_eq;
            }
        }
    }
    return a[l];
}

/// Deterministic Sn Inner Selector (O(log n))
/// For a given i, calculates the high median of |x[i] - x[j]| for j in [0, n-1], j != i.
export fn zig_sn_deterministic(x_ptr: [*]const f64, n: usize, i: usize) f64 {
    const h = n / 2; // Correct rank for j != i to match (n/2 + 1) with j = i
    return zig_sn_index_select(x_ptr, n, i, h);
}

fn zig_sn_index_select(x: [*]const f64, n: usize, i: usize, h: usize) f64 {
    const val_i = x[i];
    const na: i32 = @intCast(i);
    const nb: i32 = @intCast(n - 1 - i);

    // We want the h-th smallest in A U B.
    // A = {val_i - x[i-j] | j=1..na} (sorted)
    // B = {x[i+j] - val_i | j=1..nb} (sorted)
    
    // Standard O(log(min(na, nb))) algorithm for k-th smallest element in union of two sorted arrays.
    // However, since we access x via index midA/midB, we need to be careful with the mapping.
    
    var k = h;
    var a_start: i32 = 1;
    const a_end: i32 = na;
    var b_start: i32 = 1;
    const b_end: i32 = nb;

    while (true) {
        if (a_start > a_end) return x[i + @as(usize, @intCast(b_start + @as(i32, @intCast(k)) - 1))] - val_i;
        if (b_start > b_end) return val_i - x[i - @as(usize, @intCast(a_start + @as(i32, @intCast(k)) - 1))];
        if (k == 1) {
            const valA = val_i - x[i - @as(usize, @intCast(a_start))];
            const valB = x[i + @as(usize, @intCast(b_start))] - val_i;
            return @min(valA, valB);
        }

        const half = k / 2;
        const stepA = @min(@as(i32, @intCast(half)), a_end - a_start + 1);
        const stepB = @min(@as(i32, @intCast(half)), b_end - b_start + 1);
        
        const valA = val_i - x[i - @as(usize, @intCast(a_start + stepA - 1))];
        const valB = x[i + @as(usize, @intCast(b_start + stepB - 1))] - val_i;
        
        if (valA < valB) {
            k -= @intCast(stepA);
            a_start += stepA;
        } else {
            k -= @intCast(stepB);
            b_start += stepB;
        }
    }
}

/// Advanced Qn Selector (Johnson-Mizoguchi)
/// Returns the k-th smallest difference x[i] - x[j] for j < i.
export fn zig_qn_jm_select(x_ptr: [*]const f64, n: usize, k: u64, work_ptr: [*]f64, iweight_ptr: [*]i32, left_ptr: [*]i32, right_ptr: [*]i32) f64 {
    const x = x_ptr[0..n];
    const left = left_ptr[0..n];
    const right = right_ptr[0..n];
    const work = work_ptr[0..n];
    const iweight = iweight_ptr[0..n];

    for (0..n) |idx| {
        left[idx] = 1;
        right[idx] = @intCast(idx);
    }

    var nL: u64 = 0;
    var nR: u64 = @intCast(n * (n - 1) / 2);

    while (nR - nL > n) {
        var m: usize = 0;
        for (1..n) |idx| {
            if (left[idx] <= right[idx]) {
                const w = right[idx] - left[idx] + 1;
                const jj = left[idx] + @divFloor(w, 2);
                work[m] = x[idx] - x[idx - @as(usize, @intCast(jj))];
                iweight[m] = w;
                m += 1;
            }
        }

        const trial = whimed(work[0..m], iweight[0..m], m, @intCast(@divFloor(nR - nL, 2)));

        var sumP: u64 = 0;
        var sumQ: u64 = 0;
        
        var jp_up: usize = 0;
        var jq_up: usize = 0;
        for (1..n) |idx| {
            while (jp_up < idx and (x[idx] - x[jp_up]) >= trial) jp_up += 1;
            sumP += (idx - jp_up);
            while (jq_up < idx and (x[idx] - x[jq_up]) > trial) jq_up += 1;
            sumQ += (idx - jq_up);
        }

        if (k <= sumP) {
            var jp_refine: usize = 0;
            for (1..n) |idx| {
                while (jp_refine < idx and (x[idx] - x[jp_refine]) >= trial) jp_refine += 1;
                const jj_bound: i32 = @intCast(idx - jp_refine);
                if (jj_bound < right[idx]) {
                    right[idx] = jj_bound;
                }
            }
            nR = sumP;
        } else if (k > sumQ) {
            var jq_refine: usize = 0;
            for (1..n) |idx| {
                while (jq_refine < idx and (x[idx] - x[jq_refine]) > trial) jq_refine += 1;
                const jj_bound: i32 = @intCast(idx - jq_refine + 1);
                if (jj_bound > left[idx]) {
                    left[idx] = jj_bound;
                }
            }
            nL = sumQ;
        } else {
            return trial;
        }
    }

    var m_final: usize = 0;
    for (1..n) |idx| {
        var jj = left[idx];
        while (jj <= right[idx]) : (jj += 1) {
            work[m_final] = x[idx] - x[idx - @as(usize, @intCast(jj))];
            m_final += 1;
        }
    }

    const final_slice = work[0..m_final];
    std.sort.pdq(f64, final_slice, {}, std.sort.asc(f64));
    return final_slice[k - nL - 1];
}
