// ---------------------------------------------------------------------------
// Geometry: RDGeom::Point2D and RDGeom::Transform2D

#[derive(Clone, Copy, Debug, Default, PartialEq)]
struct P {
    x: f64,
    y: f64,
}

impl P {
    const fn new(x: f64, y: f64) -> P {
        P { x, y }
    }
    fn add(self, o: P) -> P {
        P::new(self.x + o.x, self.y + o.y)
    }
    fn sub(self, o: P) -> P {
        P::new(self.x - o.x, self.y - o.y)
    }
    fn scale(self, s: f64) -> P {
        P::new(self.x * s, self.y * s)
    }
    /// `operator-()`: multiplies by -1.0.
    fn neg(self) -> P {
        P::new(-self.x, -self.y)
    }
    fn length_sq(self) -> f64 {
        self.x * self.x + self.y * self.y
    }
    fn length(self) -> f64 {
        (self.x * self.x + self.y * self.y).sqrt()
    }
    fn dot(self, o: P) -> f64 {
        self.x * o.x + self.y * o.y
    }
    fn normalize(&mut self) -> DResult<()> {
        let ln = self.length();
        if ln < ZERO_TOLERANCE {
            return derr("Cannot normalize a zero length vector");
        }
        self.x /= ln;
        self.y /= ln;
        Ok(())
    }
}

/// `RDGeom::Transform2D`: a row-major 3x3 matrix.
#[derive(Clone, Copy, Debug)]
struct T2 {
    d: [f64; 9],
}

impl T2 {
    fn identity() -> T2 {
        let mut d = [0.0; 9];
        d[0] = 1.0;
        d[4] = 1.0;
        d[8] = 1.0;
        T2 { d }
    }

    fn set_translation(&mut self, pt: P) {
        self.d[2] = pt.x;
        self.d[5] = pt.y;
        self.d[8] = 1.0;
    }

    /// `SquareMatrix::operator*=`: `self = self * b`.
    fn mul_assign(&mut self, b: &T2) {
        let mut nd = [0.0f64; 9];
        for i in 0..3 {
            for j in 0..3 {
                let mut acc = 0.0f64;
                for k in 0..3 {
                    acc += self.d[i * 3 + k] * b.d[k * 3 + j];
                }
                nd[i * 3 + j] = acc;
            }
        }
        self.d = nd;
    }

    fn apply(&self, pt: &mut P) {
        let d = &self.d;
        let x = d[0] * pt.x + d[1] * pt.y + d[2];
        let y = d[3] * pt.x + d[4] * pt.y + d[5];
        pt.x = x;
        pt.y = y;
    }

    /// `SetTransform(pt, angle)`: rotation by `angle` about `pt`.
    fn rotation_about(pt: P, angle: f64) -> T2 {
        let mut this = T2::identity();
        let mut trans1 = T2::identity();
        trans1.set_translation(pt.neg());
        this.d[0] = angle.cos();
        this.d[1] = -angle.sin();
        this.d[3] = angle.sin();
        this.d[4] = angle.cos();
        this.mul_assign(&trans1);
        let mut trans2 = T2::identity();
        trans2.set_translation(pt);
        trans2.mul_assign(&this);
        trans2
    }

    /// `SetTransform(ref1, ref2, pt1, pt2)`.
    fn align(ref1: P, ref2: P, pt1: P, pt2: P) -> T2 {
        let rvec = ref2.sub(ref1);
        let pvec = pt2.sub(pt1);
        let dp = rvec.dot(pvec);
        let lp = rvec.length() * pvec.length();
        if lp <= 0.0 {
            return T2::identity();
        }
        let cval = (dp / lp).clamp(-1.0, 1.0);
        let mut ang = cval.acos();
        let cross = pvec.x * rvec.y - pvec.y * rvec.x;
        if cross < 0.0 {
            ang *= -1.0;
        }
        let mut this = T2::identity();
        this.d[0] = ang.cos();
        this.d[1] = -ang.sin();
        this.d[3] = ang.sin();
        this.d[4] = ang.cos();
        let mut npt1 = pt1;
        this.apply(&mut npt1);
        this.d[2] = ref1.x - npt1.x;
        this.d[5] = ref1.y - npt1.y;
        this
    }
}

/// `reflectPoint`.
fn reflect_point(point: P, loc1: P, loc2: P) -> P {
    let org = P::new(0.0, 0.0);
    let xaxis = P::new(1.0, 0.0);
    let cent = loc1.add(loc2).scale(0.5);
    let trans = T2::align(org, xaxis, cent, loc1);
    let itrans = T2::align(cent, loc1, org, xaxis);
    let mut res = point;
    trans.apply(&mut res);
    res.y = -res.y;
    itrans.apply(&mut res);
    res
}

/// `computeBisectPoint`.
fn compute_bisect_point(rcr: P, ang: f64, nb1: P, nb2: P) -> P {
    let mut cloc = nb1;
    cloc = cloc.add(nb2);
    cloc = cloc.scale(0.5);
    if ang > PI {
        cloc = cloc.sub(rcr);
        cloc = cloc.scale(-1.0);
        cloc = cloc.add(rcr);
    }
    cloc
}

/// `rotationDir`.
fn rotation_dir(center: P, loc1: P, loc2: P, rem_angle: f64) -> i32 {
    let pt1 = loc1.sub(center);
    let pt2 = loc2.sub(center);
    let mut cross = pt1.x * pt2.y - pt1.y * pt2.x;
    let diff_angle = PI - rem_angle;
    cross *= diff_angle;
    if cross >= 0.0 { -1 } else { 1 }
}

/// `computeNormal(center, other)`.
fn compute_normal(center: P, other: P) -> DResult<P> {
    let mut res = other.sub(center);
    res.normalize()?;
    Ok(P::new(-res.y, res.x))
}

/// `Point2D::angleTo`.
fn angle_to_2d(a: P, b: P) -> DResult<f64> {
    let mut t1 = a;
    let mut t2 = b;
    t1.normalize()?;
    t2.normalize()?;
    let d = t1.dot(t2).clamp(-1.0, 1.0);
    Ok(d.acos())
}

/// `computeAngle(center, loc1, loc2)`.
fn compute_angle(center: P, loc1: P, loc2: P) -> DResult<f64> {
    angle_to_2d(loc1.sub(center), loc2.sub(center))
}

/// `computeSubAngle`.
fn compute_sub_angle(degree: usize, htype: Hybridization) -> f64 {
    match htype {
        Hybridization::Unspecified | Hybridization::Sp3 => {
            if degree == 4 {
                PI / 2.0
            } else {
                2.0 * PI / 3.0
            }
        }
        Hybridization::Sp2 => 2.0 * PI / 3.0,
        _ => 2.0 * PI / degree as f64,
    }
}

/// `embedRing`.
fn embed_ring(ring: &[usize]) -> BTreeMap<usize, P> {
    let na = ring.len();
    let ang = 2.0 * PI / na as f64;
    let al = BOND_LEN / (2.0 * (1.0 - ang.cos())).sqrt();
    let mut res = BTreeMap::new();
    for (i, &a) in ring.iter().enumerate() {
        let x = al * (i as f64 * ang).cos();
        let y = al * (i as f64 * ang).sin();
        res.insert(a, P::new(x, y));
    }
    res
}
