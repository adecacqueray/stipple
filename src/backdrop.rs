//! The photo as a pale backdrop behind the marks, coloured by the same sliders as the marks.

/// Recolour `src` (RGBA) into `dst` (RGBA, same length):
/// - `tint` blends each pixel from its grey (0) to its own colour (1), like the Colour slider;
/// - `gain` multiplies the result and clips at 255, like the Contrast slider;
/// - the result is then faded towards the paper colour, keeping `strength` of it (0 = blank paper).
pub fn recolour(src: &[u8], dst: &mut [u8], tint: f32, gain: f32, strength: f32, paper: u8) {
    let p = paper as f32;
    for (s, d) in src.chunks_exact(4).zip(dst.chunks_exact_mut(4)) {
        let grey = crate::image::luma(s[0], s[1], s[2]) as f32;
        for c in 0..3 {
            let v = ((grey + (s[c] as f32 - grey) * tint) * gain).min(255.0);
            d[c] = (p + (v - p) * strength).round().clamp(0.0, 255.0) as u8;
        }
        d[3] = 255;
    }
}

#[cfg(test)]
mod tests {
    use super::recolour;

    #[test]
    fn strength_zero_is_paper_and_one_is_photo() {
        let src = [200, 40, 10, 255];
        let mut d = [0u8; 4];
        recolour(&src, &mut d, 1.0, 1.0, 0.0, 255);
        assert_eq!(d, [255, 255, 255, 255]);
        recolour(&src, &mut d, 1.0, 1.0, 1.0, 255);
        assert_eq!(d, src);
    }

    #[test]
    fn tint_zero_is_grey_and_gain_clips() {
        let mut d = [0u8; 4];
        recolour(&[200, 40, 10, 255], &mut d, 0.0, 1.0, 1.0, 0);
        assert!(d[0] == d[1] && d[1] == d[2]);
        recolour(&[200, 40, 10, 255], &mut d, 1.0, 3.0, 1.0, 0);
        assert_eq!(d[0], 255);
    }
}
