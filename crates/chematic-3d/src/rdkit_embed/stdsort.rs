//! libstdc++'s `std::sort` (introsort + final insertion sort), so that
//! unstable sorts ported from RDKit order equal elements as RDKit does.

const THRESHOLD: usize = 16;

fn lg(n: usize) -> usize {
    (usize::BITS - 1 - n.leading_zeros()) as usize
}

/// `std::sort(v.begin(), v.end(), less)`.
pub(crate) fn std_sort<T: Clone>(v: &mut [T], less: &impl Fn(&T, &T) -> bool) {
    let n = v.len();
    if n == 0 {
        return;
    }
    introsort_loop(v, 0, n, 2 * lg(n), less);
    final_insertion_sort(v, 0, n, less);
}

fn introsort_loop<T: Clone>(
    v: &mut [T],
    first: usize,
    mut last: usize,
    mut depth: usize,
    less: &impl Fn(&T, &T) -> bool,
) {
    while last - first > THRESHOLD {
        if depth == 0 {
            heap_sort(&mut v[first..last], less);
            return;
        }
        depth -= 1;
        let cut = unguarded_partition_pivot(v, first, last, less);
        introsort_loop(v, cut, last, depth, less);
        last = cut;
    }
}

fn move_median_to_first<T>(
    v: &mut [T],
    result: usize,
    a: usize,
    b: usize,
    c: usize,
    less: &impl Fn(&T, &T) -> bool,
) {
    if less(&v[a], &v[b]) {
        if less(&v[b], &v[c]) {
            v.swap(result, b);
        } else if less(&v[a], &v[c]) {
            v.swap(result, c);
        } else {
            v.swap(result, a);
        }
    } else if less(&v[a], &v[c]) {
        v.swap(result, a);
    } else if less(&v[b], &v[c]) {
        v.swap(result, c);
    } else {
        v.swap(result, b);
    }
}

fn unguarded_partition_pivot<T>(
    v: &mut [T],
    first: usize,
    last: usize,
    less: &impl Fn(&T, &T) -> bool,
) -> usize {
    let mid = first + (last - first) / 2;
    move_median_to_first(v, first, first + 1, mid, last - 1, less);
    // __unguarded_partition(first + 1, last, first)
    let pivot = first;
    let mut f = first + 1;
    let mut l = last;
    loop {
        while less(&v[f], &v[pivot]) {
            f += 1;
        }
        l -= 1;
        while less(&v[pivot], &v[l]) {
            l -= 1;
        }
        if f >= l {
            return f;
        }
        v.swap(f, l);
        f += 1;
    }
}

fn unguarded_linear_insert<T: Clone>(v: &mut [T], mut last: usize, less: &impl Fn(&T, &T) -> bool) {
    let val = v[last].clone();
    let mut next = last - 1;
    while less(&val, &v[next]) {
        v[last] = v[next].clone();
        last = next;
        next -= 1;
    }
    v[last] = val;
}

fn insertion_sort<T: Clone>(
    v: &mut [T],
    first: usize,
    last: usize,
    less: &impl Fn(&T, &T) -> bool,
) {
    if first == last {
        return;
    }
    for i in first + 1..last {
        if less(&v[i], &v[first]) {
            let val = v[i].clone();
            for k in (first..i).rev() {
                v[k + 1] = v[k].clone();
            }
            v[first] = val;
        } else {
            unguarded_linear_insert(v, i, less);
        }
    }
}

fn final_insertion_sort<T: Clone>(
    v: &mut [T],
    first: usize,
    last: usize,
    less: &impl Fn(&T, &T) -> bool,
) {
    if last - first > THRESHOLD {
        insertion_sort(v, first, first + THRESHOLD, less);
        for i in first + THRESHOLD..last {
            unguarded_linear_insert(v, i, less);
        }
    } else {
        insertion_sort(v, first, last, less);
    }
}

/// `std::__partial_sort(first, last, last)`: `make_heap` then `sort_heap`.
fn heap_sort<T: Clone>(v: &mut [T], less: &impl Fn(&T, &T) -> bool) {
    let len = v.len();
    if len < 2 {
        return;
    }
    // __make_heap
    let mut parent = (len - 2) / 2;
    loop {
        let value = v[parent].clone();
        adjust_heap(v, parent, len, value, less);
        if parent == 0 {
            break;
        }
        parent -= 1;
    }
    // __sort_heap
    let mut last = len;
    while last > 1 {
        last -= 1;
        let value = v[last].clone();
        v[last] = v[0].clone();
        adjust_heap(v, 0, last, value, less);
    }
}

fn adjust_heap<T: Clone>(
    v: &mut [T],
    mut hole: usize,
    len: usize,
    value: T,
    less: &impl Fn(&T, &T) -> bool,
) {
    let top = hole;
    let mut second = hole;
    while second < (len.saturating_sub(1)) / 2 {
        second = 2 * (second + 1);
        if less(&v[second], &v[second - 1]) {
            second -= 1;
        }
        v[hole] = v[second].clone();
        hole = second;
    }
    if len & 1 == 0 && second == (len - 2) / 2 {
        second = 2 * (second + 1);
        v[hole] = v[second - 1].clone();
        hole = second - 1;
    }
    // __push_heap
    let mut parent = if hole > 0 { (hole - 1) / 2 } else { 0 };
    while hole > top && less(&v[parent], &value) {
        v[hole] = v[parent].clone();
        hole = parent;
        parent = if hole > 0 { (hole - 1) / 2 } else { 0 };
    }
    v[hole] = value;
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sorts_like_a_sort() {
        let mut x: Vec<(usize, usize)> = (0..200).map(|i| ((i * 7919) % 13, i)).collect();
        std_sort(&mut x, &|a, b| a.0 < b.0);
        assert!(x.windows(2).all(|w| w[0].0 <= w[1].0));
        let mut y: Vec<(usize, usize)> = (0..10).map(|i| (i % 3, i)).collect();
        std_sort(&mut y, &|a, b| a.0 < b.0);
        // insertion sort (n <= 16) is stable
        assert_eq!(
            y,
            vec![
                (0, 0),
                (0, 3),
                (0, 6),
                (0, 9),
                (1, 1),
                (1, 4),
                (1, 7),
                (2, 2),
                (2, 5),
                (2, 8)
            ]
        );
    }
}
