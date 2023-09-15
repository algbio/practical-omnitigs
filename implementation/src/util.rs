use log::debug;

pub fn median<T: Ord>(slice: &mut [T]) -> Option<&mut T> {
    debug!("Computing median of {} elements", slice.len());

    if slice.is_empty() {
        None
    } else {
        Some(slice.select_nth_unstable(slice.len() / 2).1)
    }
}

pub fn mean(slice: &[f64]) -> Option<f64> {
    debug!("Computing mean of {} elements", slice.len());

    if slice.is_empty() {
        None
    } else {
        Some(slice.iter().cloned().sum::<f64>() / slice.len() as f64)
    }
}
