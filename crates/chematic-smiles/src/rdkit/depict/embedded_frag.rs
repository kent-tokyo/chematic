// ---------------------------------------------------------------------------
// EmbeddedAtom / EmbeddedFrag

#[derive(Clone, Debug)]
struct EAtom {
    angle: f64,
    nbr1: i64,
    nbr2: i64,
    cis_trans_nbr: i64,
    ccw: bool,
    rot_dir: i32,
    loc: P,
    normal: P,
    neighs: Vec<usize>,
    density: f64,
    fixed: bool,
}

impl Default for EAtom {
    fn default() -> EAtom {
        EAtom {
            angle: -1.0,
            nbr1: -1,
            nbr2: -1,
            cis_trans_nbr: -1,
            ccw: true,
            rot_dir: 0,
            loc: P::default(),
            normal: P::default(),
            neighs: Vec::new(),
            density: -1.0,
            fixed: false,
        }
    }
}

impl EAtom {
    fn transform(&mut self, t: &T2) {
        let mut temp = self.loc.add(self.normal);
        t.apply(&mut self.loc);
        t.apply(&mut temp);
        self.normal = temp.sub(self.loc);
    }
}

fn key(i: i64) -> DResult<usize> {
    if i < 0 {
        derr("negative atom index used as a map key")
    } else {
        Ok(i as usize)
    }
}

#[derive(Clone, Debug, Default)]
struct EFrag {
    done: bool,
    /// Removed from the `std::list` of fragments.
    dead: bool,
    px: f64,
    nx: f64,
    py: f64,
    ny: f64,
    ea: BTreeMap<usize, EAtom>,
    /// `d_attachPts` (a `std::list`).
    attach: Vec<usize>,
}

impl EFrag {
    /// `d_eatoms[aid]` (inserting a default atom when absent).
    fn at(&mut self, aid: usize) -> &mut EAtom {
        self.ea.entry(aid).or_default()
    }

    /// `d_eatoms.at(aid)`.
    fn get(&self, aid: usize) -> DResult<&EAtom> {
        match self.ea.get(&aid) {
            Some(a) => Ok(a),
            None => derr("embedded atom not found"),
        }
    }

    fn loc(&self, aid: usize) -> DResult<P> {
        Ok(self.get(aid)?.loc)
    }

    fn has(&self, aid: usize) -> bool {
        self.ea.contains_key(&aid)
    }

    /// `EmbeddedFrag(aid, mol)`.
    fn from_atom(ctx: &Ctx, aid: usize) -> DResult<EFrag> {
        let mut f = EFrag::default();
        let eatm = EAtom {
            loc: P::new(0.0, 0.0),
            normal: P::new(1.0, 0.0),
            angle: -1.0,
            ccw: true,
            ..EAtom::default()
        };
        f.ea.insert(aid, eatm);
        f.update_new_neighs(ctx, aid)?;
        Ok(f)
    }

    /// `EmbeddedFrag(mol, coordMap)`: fixed atoms at the given positions.
    fn from_coord_map(ctx: &Ctx, map: &BTreeMap<usize, P>) -> DResult<EFrag> {
        let mut f = EFrag::default();
        for (&aid, &loc) in map {
            let eatm = EAtom {
                loc,
                fixed: true,
                ..EAtom::default()
            };
            f.ea.insert(aid, eatm);
        }
        f.setup_new_neighs(ctx)?;
        f.setup_attachment_points(ctx)?;
        Ok(f)
    }

    /// `setupAttachmentPoints`.
    fn setup_attachment_points(&mut self, ctx: &Ctx) -> DResult<()> {
        for dai in self.attach.clone() {
            let enbrs = self.get(dai)?.neighs.clone();
            let done: Vec<usize> = ctx.mol.nbrs(dai).filter(|nb| !enbrs.contains(nb)).collect();
            match done.len() {
                0 => {
                    let a = self.at(dai);
                    a.normal = P::new(1.0, 0.0);
                    a.angle = -1.0;
                }
                1 => {
                    let nbid = done[0];
                    let normal = compute_normal(self.loc(dai)?, self.loc(nbid)?)?;
                    let a = self.at(dai);
                    a.nbr1 = nbid as i64;
                    a.normal = normal;
                }
                2 => {
                    let ang =
                        compute_angle(self.loc(dai)?, self.loc(done[0])?, self.loc(done[1])?)?;
                    let a = self.at(dai);
                    a.nbr1 = done[0] as i64;
                    a.nbr2 = done[1] as i64;
                    a.angle = ang;
                }
                _ => self.compute_nbrs_and_ang(ctx, dai, &done)?,
            }
        }
        Ok(())
    }

    /// `computeNbrsAndAng`.
    fn compute_nbrs_and_ang(&mut self, ctx: &Ctx, aid: usize, done: &[usize]) -> DResult<()> {
        let center = self.loc(aid)?;
        let mut pairs: Vec<(f64, (usize, usize))> = Vec::new();
        // RDKit's inner loop starts at `nbi3++` (a post-increment), so it
        // pairs every neighbour with itself too.
        for i in 0..done.len() {
            for j in i..done.len() {
                let ang = compute_angle(center, self.loc(done[i])?, self.loc(done[j])?)?;
                pairs.push((ang, (done[i], done[j])));
            }
        }
        // `std::list::sort` is stable.
        pairs.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
        let mut winner = *pairs.last().expect("three neighbours");
        let ri = ctx.mol.ring_info();
        for pr in pairs.iter().rev() {
            if ri.num_atom_rings(pr.1.0) <= 1 && ri.num_atom_rings(pr.1.1) <= 1 {
                winner = *pr;
                break;
            }
        }
        let (wnb1, wnb2) = winner.1;
        let (mut nb1, mut nb2) = (-1i64, -1i64);
        for &(_, (f, s)) in &pairs {
            if wnb1 == f {
                nb2 = wnb1 as i64;
                nb1 = s as i64;
                break;
            } else if wnb1 == s {
                nb2 = wnb1 as i64;
                nb1 = f as i64;
                break;
            } else if wnb2 == f {
                nb2 = wnb2 as i64;
                nb1 = s as i64;
                break;
            } else if wnb2 == s {
                nb2 = wnb2 as i64;
                nb1 = f as i64;
                break;
            }
        }
        let w_ang = winner.0;
        let rot = rotation_dir(center, self.loc(key(nb1)?)?, self.loc(key(nb2)?)?, w_ang);
        let a = self.at(aid);
        a.rot_dir = rot;
        a.nbr1 = nb1;
        a.nbr2 = nb2;
        a.angle = 2.0 * PI - w_ang;
        Ok(())
    }

    /// `EmbeddedFrag(dblBond)`.
    fn from_double_bond(b: usize, ctx: &Ctx) -> EFrag {
        let bond = &ctx.mol.bonds[b];
        let beg = bond.begin;
        let end = bond.end;
        let mut f = EFrag::default();
        let beatm = EAtom {
            loc: P::new(0.0, 0.0),
            nbr1: end as i64,
            normal: P::new(0.0, -1.0),
            ccw: false,
            cis_trans_nbr: bond.stereo_atoms[0] as i64,
            ..EAtom::default()
        };
        *f.at(beg) = beatm;
        let (normal, ccw) = if bond.stereo == BondStereo::Z {
            (P::new(0.0, -1.0), true)
        } else {
            (P::new(0.0, 1.0), false)
        };
        let eeatm = EAtom {
            loc: P::new(BOND_LEN, 0.0),
            nbr1: beg as i64,
            cis_trans_nbr: bond.stereo_atoms[1] as i64,
            normal,
            ccw,
            ..EAtom::default()
        };
        *f.at(end) = eeatm;
        f
    }

    /// `EmbeddedFrag(mol, fusedRings, false)`.
    fn from_fused_rings(ctx: &Ctx, fused: &[Vec<usize>]) -> DResult<EFrag> {
        let mut f = EFrag::default();
        f.embed_fused_rings(ctx, fused)?;
        Ok(f)
    }

    fn size(&self) -> usize {
        self.ea.len()
    }

    /// `findNumNeigh`.
    fn find_num_neigh(&self, pt: P, radius: f64) -> i32 {
        let mut res = 0;
        for a in self.ea.values() {
            if a.loc.sub(pt).length() < radius {
                res += 1;
            }
        }
        res
    }

    /// `updateNewNeighs`.
    fn update_new_neighs(&mut self, ctx: &Ctx, aid: usize) -> DResult<()> {
        let mut neighs = Vec::new();
        let mut h_indices = Vec::new();
        for nb in ctx.mol.nbrs(aid) {
            if !self.has(nb) {
                if ctx.anum(nb) != 1 {
                    neighs.push(nb);
                } else {
                    h_indices.push(nb);
                }
            }
        }
        neighs.extend(h_indices);
        let deg = ctx.degree(aid);
        if !neighs.is_empty() && (deg < 4 || neighs.len() < 3) {
            neighs = ctx.rank_atoms(&neighs);
        } else if deg >= 4 && neighs.len() >= 3 {
            neighs = ctx.set_nbr_order(aid, &neighs)?;
        }
        let nonempty = !neighs.is_empty();
        self.at(aid).neighs = neighs;
        if nonempty && !self.attach.contains(&aid) {
            self.attach.push(aid);
        }
        Ok(())
    }

    /// `setupNewNeighs`.
    fn setup_new_neighs(&mut self, ctx: &Ctx) -> DResult<()> {
        self.attach.clear();
        let keys: Vec<usize> = self.ea.keys().copied().collect();
        for k in keys {
            self.update_new_neighs(ctx, k)?;
        }
        self.attach = ctx.rank_atoms(&self.attach);
        Ok(())
    }

    /// `embedFusedRings` (without ring templates).
    fn embed_fused_rings(&mut self, ctx: &Ctx, fused: &[Vec<usize>]) -> DResult<()> {
        let mut coords: Vec<BTreeMap<usize, P>> = Vec::with_capacity(fused.len());
        for ring in fused {
            let mut rc = embed_ring(ring);
            mirror_trans_ring_atoms(ctx, ring, &mut rc)?;
            coords.push(rc);
        }
        let mut done_rings: Vec<usize> = Vec::new();
        let first = pick_first_ring_to_embed(ctx, fused);
        self.init_from_ring_coords(&fused[first], &coords[first]);
        done_rings.push(first);
        let funion = union_rings(fused, None);
        while self.ea.len() < funion.len() {
            let (common, next_id) = find_next_ring_to_embed(&done_rings, fused)?;
            let mut emb_ring = EFrag::default();
            emb_ring.init_from_ring_coords(&fused[next_id], &coords[next_id]);
            let mut pin = Vec::new();
            if common.len() == 1 {
                let t = self.compute_one_atom_trans(common[0], &emb_ring)?;
                emb_ring.transform(&t);
                pin.push(common[0]);
            } else {
                let aid1 = common[0];
                let aid2 = *common.last().unwrap();
                pin.push(aid1);
                pin.push(aid2);
                let loc1 = coords[next_id][&aid1];
                let loc2 = coords[next_id][&aid2];
                let t = T2::align(self.loc(aid1)?, self.loc(aid2)?, loc1, loc2);
                emb_ring.transform(&t);
                self.reflect_if_necessary_density(&mut emb_ring, aid1, aid2);
            }
            self.merge_ring(&emb_ring, common.len(), &pin);
            done_rings.push(next_id);
        }
        Ok(())
    }

    /// `initFromRingCoords`.
    fn init_from_ring_coords(&mut self, ring: &[usize], map: &BTreeMap<usize, P>) {
        let largest_angle = PI * (1.0 - (2.0 / ring.len() as f64));
        let mut prev = *ring.last().unwrap();
        for (cnt, &ai) in ring.iter().enumerate() {
            let eatm = EAtom {
                loc: map[&ai],
                angle: largest_angle,
                nbr1: prev as i64,
                ..EAtom::default()
            };
            if cnt > 0 {
                self.at(prev).nbr2 = ai as i64;
            }
            assign_eatom(self.at(ai), &eatm);
            prev = ai;
        }
        self.at(prev).nbr2 = ring[0] as i64;
    }

    /// `computeOneAtomTrans`.
    fn compute_one_atom_trans(&mut self, comm: usize, other: &EFrag) -> DResult<T2> {
        let rcr = self.at(comm).loc;
        let oeatm = other.get(comm)?;
        let ccr = oeatm.loc;
        let (onb1, onb2) = (oeatm.nbr1, oeatm.nbr2);
        if onb1 < 0 || onb2 < 0 {
            return derr("computeOneAtomTrans: missing neighbors");
        }
        let mut mid = other.loc(onb1 as usize)?;
        mid = mid.add(other.loc(onb2 as usize)?);
        mid = mid.scale(0.5);
        let nb1 = key(self.at(comm).nbr1)?;
        let nb2 = key(self.at(comm).nbr2)?;
        let nbp1 = self.at(nb1).loc;
        let nbp2 = self.at(nb2).loc;
        let ang = self.at(comm).angle;
        let largest = 2.0 * PI - ang;
        let bpt = compute_bisect_point(rcr, largest, nbp1, nbp2);
        Ok(T2::align(rcr, bpt, ccr, mid))
    }

    fn transform(&mut self, t: &T2) {
        for a in self.ea.values_mut() {
            a.transform(t);
        }
    }

    /// `EmbeddedFrag::Reflect(loc1, loc2)` with fixed line points.
    fn reflect(&mut self, loc1: P, loc2: P) {
        for a in self.ea.values_mut() {
            let temp = a.loc.add(a.normal);
            a.loc = reflect_point(a.loc, loc1, loc2);
            let temp = reflect_point(temp, loc1, loc2);
            a.normal = temp.sub(a.loc);
            a.ccw = !a.ccw;
        }
    }

    /// `EmbeddedAtom::Reflect(loc1, loc2)` on atom `k` where the line runs
    /// through this fragment's atoms `l1` and `l2` (the C++ holds
    /// references to their locations, which may be `k`'s own).
    fn reflect_atom_by_refs(&mut self, k: usize, l1: usize, l2: usize) -> DResult<()> {
        let loc = self.loc(k)?;
        let normal = self.get(k)?.normal;
        let temp = loc.add(normal);
        let new_loc = reflect_point(loc, self.loc(l1)?, self.loc(l2)?);
        self.at(k).loc = new_loc;
        let temp = reflect_point(temp, self.loc(l1)?, self.loc(l2)?);
        let a = self.at(k);
        a.normal = temp.sub(a.loc);
        a.ccw = !a.ccw;
        Ok(())
    }

    /// `reflectIfNecessaryDensity`.
    fn reflect_if_necessary_density(&mut self, emb: &mut EFrag, aid1: usize, aid2: usize) {
        let pin1 = self.at(aid1).loc;
        let pin2 = self.at(aid2).loc;
        let mut dn = 0.0f64;
        let mut dr = 0.0f64;
        for (k, oa) in &emb.ea {
            if !self.ea.contains_key(k) {
                let loc1 = oa.loc;
                let rloc1 = reflect_point(loc1, pin1, pin2);
                for ta in self.ea.values() {
                    let td = ta.loc.sub(loc1).length();
                    let rtd = ta.loc.sub(rloc1).length();
                    if td > 1.0e-3 {
                        dn += 1.0 / td;
                    } else {
                        dn += 1000.0;
                    }
                    if rtd > 1.0e-3 {
                        dr += 1.0 / rtd;
                    } else {
                        dr += 1000.0;
                    }
                }
            }
        }
        if dn - dr > 1.0e-4 {
            emb.reflect(pin1, pin2);
        }
    }

    /// `reflectIfNecessaryCisTrans`.
    fn reflect_if_necessary_cis_trans(
        &mut self,
        emb: &mut EFrag,
        ct_case: u32,
        aid1: usize,
        aid2: usize,
    ) -> DResult<()> {
        let p1 = self.at(aid1).loc;
        let mut ratm;
        let p1norm;
        if ct_case == 1 {
            p1norm = emb.at(aid1).normal;
            let ring_atm = emb.at(aid1).cis_trans_nbr;
            if ring_atm >= 0 && self.has(ring_atm as usize) {
                ratm = self.at(ring_atm as usize).loc;
            } else {
                // RDKit warns that the double bond stereo may be wrong.
                return Ok(());
            }
        } else {
            p1norm = self.at(aid1).normal;
            let ring_atm = key(self.at(aid1).cis_trans_nbr)?;
            ratm = emb.at(ring_atm).loc;
        }
        ratm = ratm.sub(p1);
        let dot = ratm.dot(p1norm);
        let p2 = self.at(aid2).loc;
        if dot < 0.0 {
            emb.reflect(p1, p2);
        }
        Ok(())
    }

    /// `reflectIfNecessaryThirdPt`.
    fn reflect_if_necessary_third_pt(
        &mut self,
        emb: &mut EFrag,
        aid1: usize,
        aid2: usize,
        aid3: usize,
    ) -> DResult<()> {
        let pt1 = self.at(aid1).loc;
        let pt2 = self.at(aid2).loc;
        let mut normal = pt2.sub(pt1);
        // rotate90
        let t = normal.x;
        normal.x = -normal.y;
        normal.y = t;
        let oth3 = emb.loc(aid3)?.sub(pt1);
        let pt3 = self.at(aid3).loc.sub(pt1);
        let dot1 = normal.dot(pt3);
        let dot2 = normal.dot(oth3);
        if dot1 * dot2 < 0.0 {
            emb.reflect(pt1, pt2);
        }
        Ok(())
    }

    /// `mergeRing`.
    fn merge_ring(&mut self, emb: &EFrag, n_common: usize, pin: &[usize]) {
        for (&aid, oa) in &emb.ea {
            if !self.ea.contains_key(&aid) {
                assign_eatom(self.at(aid), oa);
            } else if n_common <= 2 && pin.contains(&aid) {
                let a = self.at(aid);
                a.angle += oa.angle;
                if a.nbr1 == oa.nbr1 {
                    a.nbr1 = oa.nbr2;
                } else if a.nbr1 == oa.nbr2 {
                    a.nbr1 = oa.nbr1;
                } else if a.nbr2 == oa.nbr1 {
                    a.nbr2 = oa.nbr2;
                } else if a.nbr2 == oa.nbr2 {
                    a.nbr2 = oa.nbr1;
                }
            }
        }
    }

    /// `addNonRingAtom`.
    fn add_non_ring_atom(&mut self, ctx: &Ctx, aid: usize, to: usize) -> DResult<()> {
        if self.has(aid) || !self.has(to) {
            return derr("addNonRingAtom: bad atoms");
        }
        if self.at(to).angle > 0.0 {
            self.add_atom_to_atom_with_ang(aid, to)?;
        } else {
            self.add_atom_to_atom_with_no_ang(ctx, aid, to)?;
        }
        // neighs.erase(std::remove(...)): erases the single removed slot; with
        // `aid` absent libstdc++ drops the last element.
        let neighs = &mut self.at(to).neighs;
        if let Some(p) = neighs.iter().position(|&x| x == aid) {
            neighs.remove(p);
        } else {
            neighs.pop();
        }
        self.update_new_neighs(ctx, aid)
    }

    /// `addAtomToAtomWithAng`.
    fn add_atom_to_atom_with_ang(&mut self, aid: usize, to: usize) -> DResult<()> {
        let r = self.get(to)?.clone();
        let ref_loc = r.loc;
        let nnbr = r.neighs.len();
        let rem_angle = 2.0 * PI - r.angle;
        let mut curr_angle = rem_angle / (1 + nnbr) as f64;
        self.at(to).angle += curr_angle;
        let nb1 = self.loc(key(r.nbr1)?)?;
        let nb2 = self.loc(key(r.nbr2)?)?;
        if self.at(to).rot_dir == 0 {
            self.at(to).rot_dir = rotation_dir(ref_loc, nb1, nb2, rem_angle);
        }
        curr_angle *= f64::from(self.at(to).rot_dir);
        let rtrans = T2::rotation_about(ref_loc, curr_angle);
        let mut curr_loc = nb2;
        rtrans.apply(&mut curr_loc);
        if rem_angle.abs() - PI < 1e-3 {
            let mut curr_loc2 = nb2;
            let rtrans2 = T2::rotation_about(ref_loc, -curr_angle);
            rtrans2.apply(&mut curr_loc2);
            if self.find_num_neigh(curr_loc, 0.5) > self.find_num_neigh(curr_loc2, 0.5) {
                curr_loc = curr_loc2;
            }
        }
        self.at(to).nbr2 = aid as i64;
        let mut eatm = EAtom {
            loc: curr_loc,
            nbr1: to as i64,
            angle: -1.0,
            ..EAtom::default()
        };
        let tpt = curr_loc.sub(ref_loc);
        let mut norm = P::new(-tpt.y, tpt.x);
        let tp1 = curr_loc.add(norm);
        let tp2 = curr_loc.sub(norm);
        let nccw = self.find_num_neigh(tp1, NEIGH_RADIUS);
        let ncw = self.find_num_neigh(tp2, NEIGH_RADIUS);
        norm.normalize()?;
        if nccw < ncw {
            eatm.normal = norm;
            eatm.ccw = false;
        } else {
            eatm.normal = norm.neg();
            eatm.ccw = true;
        }
        assign_eatom(self.at(aid), &eatm);
        Ok(())
    }

    /// `addAtomToAtomWithNoAng`.
    fn add_atom_to_atom_with_no_ang(&mut self, ctx: &Ctx, aid: usize, to: usize) -> DResult<()> {
        let r = self.get(to)?.clone();
        let ref_loc = r.loc;
        let mut ref_ccw = r.ccw;
        let mut curr_loc = r.normal;
        if r.cis_trans_nbr >= 0 && r.cis_trans_nbr as usize != aid {
            ref_ccw = !ref_ccw;
            curr_loc = curr_loc.scale(-1.0);
        }
        if curr_loc.length_sq() <= 1.0e-8 {
            return derr("addAtomToAtomWithNoAng: zero normal");
        }
        let deg = ctx.degree(to);
        let mut angle = compute_sub_angle(deg, ctx.mol.atoms[to].hybrid);
        let mut flip_norm = false;
        if self.at(to).nbr1 >= 0 {
            self.at(to).angle = angle;
            self.at(to).nbr2 = aid as i64;
        } else {
            let mut norm = self.get(to)?.normal;
            let rtrans = T2::rotation_about(P::new(0.0, 0.0), angle);
            rtrans.apply(&mut norm);
            self.at(to).normal = norm;
            self.at(to).nbr1 = aid as i64;
            flip_norm = true;
        }
        angle -= PI / 2.0;
        if !ref_ccw {
            angle *= -1.0;
        }
        let trans = T2::rotation_about(P::new(0.0, 0.0), angle);
        trans.apply(&mut curr_loc);
        curr_loc = curr_loc.scale(BOND_LEN);
        curr_loc = curr_loc.add(ref_loc);
        let tpt = ref_loc.sub(curr_loc);
        let mut norm = P::new(-tpt.y, tpt.x);
        if ref_ccw ^ flip_norm {
            norm = norm.scale(-1.0);
        }
        norm.normalize()?;
        let eatm = EAtom {
            loc: curr_loc,
            normal: norm,
            nbr1: to as i64,
            angle: -1.0,
            ccw: (!ref_ccw) ^ flip_norm,
            ..EAtom::default()
        };
        assign_eatom(self.at(aid), &eatm);
        Ok(())
    }

    /// `findCommonAtoms`.
    fn find_common_atoms(&self, other: &EFrag) -> Vec<usize> {
        self.ea
            .keys()
            .filter(|k| other.ea.contains_key(k))
            .copied()
            .collect()
    }

    /// `mergeNoCommon`.
    fn merge_no_common(
        &mut self,
        ctx: &Ctx,
        emb: &mut EFrag,
        to: usize,
        nbr: usize,
    ) -> DResult<()> {
        self.add_non_ring_atom(ctx, nbr, to)?;
        emb.add_non_ring_atom(ctx, to, nbr)?;
        let mut comm = vec![to, nbr];
        self.merge_with_common(ctx, emb, &mut comm)
    }

    /// `mergeWithCommon`.
    fn merge_with_common(
        &mut self,
        ctx: &Ctx,
        emb: &mut EFrag,
        comm: &mut Vec<usize>,
    ) -> DResult<()> {
        let mut ct_case = 0u32;
        if comm.len() == 1 {
            let c = comm[0];
            let other: i64;
            if self.at(c).cis_trans_nbr >= 0 {
                ct_case = 2;
                other = self.at(c).nbr1;
                emb.add_non_ring_atom(ctx, key(other)?, c)?;
            } else if emb.at(c).cis_trans_nbr >= 0 {
                ct_case = 1;
                other = emb.at(c).nbr1;
                self.add_non_ring_atom(ctx, key(other)?, c)?;
            } else {
                other = self.at(c).nbr1;
                if other >= 0 {
                    emb.add_non_ring_atom(ctx, other as usize, c)?;
                }
            }
            if other >= 0 {
                comm.push(other as usize);
            }
        }
        let rtrans = if comm.len() == 1 {
            self.compute_one_atom_trans(comm[0], emb)?
        } else {
            let (c1, c2) = (comm[0], comm[1]);
            T2::align(self.loc(c1)?, self.loc(c2)?, emb.loc(c1)?, emb.loc(c2)?)
        };
        emb.transform(&rtrans);
        if comm.len() >= 2 {
            if ct_case > 0 {
                self.reflect_if_necessary_cis_trans(emb, ct_case, comm[0], comm[1])?;
            } else if comm.len() == 2 {
                self.reflect_if_necessary_density(emb, comm[0], comm[1]);
            } else {
                self.reflect_if_necessary_third_pt(emb, comm[0], comm[1], comm[2])?;
            }
        }
        self.merge_embedded_atoms(emb, comm);
        for &c in comm.iter() {
            self.update_new_neighs(ctx, c)?;
        }
        Ok(())
    }

    /// `mergeFragsWithComm`.
    fn merge_frags_with_comm(&mut self, ctx: &Ctx, efrags: &mut [EFrag]) -> DResult<()> {
        loop {
            let mut comm = Vec::new();
            let mut nfri = None;
            for (i, f) in efrags.iter().enumerate() {
                if !f.dead && !f.done {
                    comm = self.find_common_atoms(f);
                    if !comm.is_empty() {
                        nfri = Some(i);
                        break;
                    }
                }
            }
            if comm.is_empty() {
                break;
            }
            let i = nfri.unwrap();
            let mut other = std::mem::take(&mut efrags[i]);
            self.merge_with_common(ctx, &mut other, &mut comm)?;
            for &c in &comm {
                if self.get(c)?.neighs.is_empty()
                    && let Some(p) = self.attach.iter().position(|&x| x == c)
                {
                    self.attach.remove(p);
                }
            }
            efrags[i].dead = true;
        }
        Ok(())
    }

    /// `expandEfrag`.
    fn expand(&mut self, ctx: &Ctx, nratms: &mut Vec<usize>, efrags: &mut [EFrag]) -> DResult<()> {
        self.merge_frags_with_comm(ctx, efrags)?;
        while !self.attach.is_empty() {
            let aid = self.attach[0];
            let nbrs = self.at(aid).neighs.clone();
            if nbrs.is_empty() {
                return derr("expandEfrag: attachment point without neighbors");
            }
            for nbri in nbrs {
                if let Some(p) = nratms.iter().position(|&x| x == nbri) {
                    self.add_non_ring_atom(ctx, nbri, aid)?;
                    nratms.remove(p);
                } else {
                    let mut nfri = None;
                    for (i, f) in efrags.iter().enumerate() {
                        if !f.dead && !f.done && f.ea.contains_key(&nbri) {
                            nfri = Some(i);
                            break;
                        }
                    }
                    if let Some(i) = nfri {
                        let mut other = std::mem::take(&mut efrags[i]);
                        self.merge_no_common(ctx, &mut other, aid, nbri)?;
                        if self.get(nbri)?.neighs.is_empty()
                            && let Some(p) = self.attach.iter().position(|&x| x == nbri)
                        {
                            self.attach.remove(p);
                        }
                        efrags[i].dead = true;
                    }
                }
            }
            self.attach.remove(0);
            self.at(aid).neighs.clear();
            self.merge_frags_with_comm(ctx, efrags)?;
        }
        Ok(())
    }

    /// `computeBox`.
    fn compute_box(&mut self) {
        self.px = -1.0e8;
        self.nx = 1.0e8;
        self.py = -1.0e8;
        self.ny = 1.0e8;
        for a in self.ea.values() {
            self.px = self.px.max(a.loc.x);
            self.nx = self.nx.min(a.loc.x);
            self.py = self.py.max(a.loc.y);
            self.ny = self.ny.min(a.loc.y);
        }
        self.nx *= -1.0;
        self.ny *= -1.0;
    }

    /// `canonicalizeOrientation`.
    fn canonicalize_orientation(&mut self) -> DResult<()> {
        if self.ea.len() <= 1 {
            return Ok(());
        }
        let mut cent = P::new(0.0, 0.0);
        for a in self.ea.values() {
            cent = cent.add(a.loc);
        }
        cent = cent.scale(1.0 / self.ea.len() as f64);
        let (mut xx, mut xy, mut yy) = (0.0f64, 0.0f64, 0.0f64);
        for a in self.ea.values_mut() {
            a.loc = a.loc.sub(cent);
            xx += a.loc.x * a.loc.x;
            xy += a.loc.x * a.loc.y;
            yy += a.loc.y * a.loc.y;
        }
        let mut d = (xx - yy) * (xx - yy) + 4.0 * xy * xy;
        d = d.sqrt();
        let mut eig1 = P::new(2.0 * xy, (yy - xx) + d);
        if eig1.length() <= 1e-4 {
            return Ok(());
        }
        let e_val1 = (xx + yy + d) / 2.0;
        eig1.normalize()?;
        let mut eig2 = P::new(2.0 * xy, (yy - xx) - d);
        let e_val2 = (xx + yy - d) / 2.0;
        if eig2.length() > 1e-4 {
            eig2.normalize()?;
            if e_val2 > e_val1 {
                std::mem::swap(&mut eig1, &mut eig2);
            }
        }
        let mut trans = T2::identity();
        trans.d[0] = eig1.x;
        trans.d[3] = -eig1.y;
        trans.d[1] = eig1.y;
        trans.d[4] = eig1.x;
        self.transform(&trans);
        Ok(())
    }

    /// `findCollisions`.
    fn find_collisions(&mut self, ctx: &Ctx, include_bonds: bool) -> Vec<(usize, usize)> {
        let mut res = Vec::new();
        for a in self.ea.values_mut() {
            a.density = 0.0;
        }
        let col_thres2 = COLLISION_THRES * COLLISION_THRES;
        let keys: Vec<usize> = self.ea.keys().copied().collect();
        let locs: Vec<P> = self.ea.values().map(|a| a.loc).collect();
        let mut dens = vec![0.0f64; keys.len()];
        for i in 0..keys.len() {
            let f1 = if ctx.anum(keys[i]) != 6 {
                HETEROATOM_COLL_SCALE
            } else {
                1.0
            };
            for j in 0..i {
                let f2 = if ctx.anum(keys[j]) != 6 {
                    HETEROATOM_COLL_SCALE
                } else {
                    1.0
                };
                let mut d2 = locs[j].sub(locs[i]).length_sq();
                if d2 > 1.0e-3 {
                    dens[i] += 1.0 / d2;
                    dens[j] += 1.0 / d2;
                } else {
                    dens[i] += 1000.0;
                    dens[j] += 1000.0;
                }
                d2 /= f1 * f2;
                if d2 < col_thres2 {
                    res.push((keys[i], keys[j]));
                }
            }
        }
        for (a, d) in self.ea.values_mut().zip(dens) {
            a.density = d;
        }
        if include_bonds {
            let bond_thres2 = BOND_THRES * BOND_THRES;
            let bonds = &ctx.mol.bonds;
            for (bid1, b1) in bonds.iter().enumerate() {
                let (beg1, end1) = (b1.begin, b1.end);
                let (Some(l_beg1), Some(l_end1)) = (self.ea.get(&beg1), self.ea.get(&end1)) else {
                    continue;
                };
                let (l_beg1, l_end1) = (l_beg1.loc, l_end1.loc);
                let v1 = l_end1.sub(l_beg1);
                let avg1 = l_end1.add(l_beg1).scale(0.5);
                for b2 in &bonds[bid1 + 1..] {
                    let (beg2, end2) = (b2.begin, b2.end);
                    let (Some(l_beg2), Some(l_end2)) = (self.ea.get(&beg2), self.ea.get(&end2))
                    else {
                        continue;
                    };
                    let (l_beg2, l_end2) = (l_beg2.loc, l_end2.loc);
                    let avg2 = l_end2.add(l_beg2).scale(0.5).sub(avg1);
                    if avg2.length_sq() < 0.5 && avg2.length_sq() < bond_thres2 {
                        let v2 = l_beg2.sub(l_beg1);
                        let v3 = l_end2.sub(l_beg1);
                        let val_prod = cross_val(v1, v2) * cross_val(v1, v3);
                        if val_prod < -1e-6 {
                            res.push(ctx.closest_pair(beg1, end1, beg2, end2));
                        }
                    }
                }
            }
        }
        res
    }

    /// `totalDensity`.
    fn total_density(&self) -> f64 {
        let mut acc = 0.0f64;
        for a in self.ea.values() {
            acc += a.density;
        }
        acc
    }

    /// `flipAboutBond`.
    fn flip_about_bond(&mut self, ctx: &Ctx, bond_id: usize, flip_end: bool) -> DResult<()> {
        if ctx.num_bond_rings(bond_id) > 0 {
            return derr("flipAboutBond: ring bond");
        }
        let bond = &ctx.mol.bonds[bond_id];
        let (mut beg, mut end) = (bond.begin, bond.end);
        if !flip_end {
            std::mem::swap(&mut beg, &mut end);
        }
        self.get(beg)?;
        self.get(end)?;
        let mut end_side = Vec::new();
        ctx.recurse_atom_one_side(end, beg, &mut end_side);
        let n_fixed = self.ea.values().filter(|a| a.fixed).count();
        let mut n_end_fixed = 0;
        if n_fixed > 0 {
            for &e in &end_side {
                if self.at(e).fixed {
                    n_end_fixed += 1;
                }
            }
        }
        let mut end_side_flip = true;
        if n_end_fixed > 0 {
            return Ok(());
        } else {
            let nats = self.ea.len();
            let n_end = end_side.len();
            if nats.wrapping_sub(n_end) < n_end {
                end_side_flip = false;
            }
        }
        let keys: Vec<usize> = self.ea.keys().copied().collect();
        for k in keys {
            let in_end = end_side.contains(&k);
            if end_side_flip ^ !in_end {
                self.reflect_atom_by_refs(k, beg, end)?;
            }
        }
        Ok(())
    }

    /// `openAngles`.
    fn open_angles(&mut self, ctx: &Ctx, aid1: usize, aid2: usize) -> DResult<()> {
        let deg1 = ctx.degree(aid1);
        let deg2 = ctx.degree(aid2);
        let fixed1 = self.get(aid1)?.fixed;
        let fixed2 = self.get(aid2)?.fixed;
        if (deg1 > 1 || fixed1) && (deg2 > 1 || fixed2) {
            return Ok(());
        }
        let (aid_a, aid_b, typ);
        if (deg1 == 1 && !fixed1) && (deg2 == 1 && !fixed2) {
            aid_a = ctx.deg1_neighbor(aid1)?;
            aid_b = ctx.deg1_neighbor(aid2)?;
            typ = 1;
        } else if (deg1 == 1 && !fixed1) && (deg2 > 1 || fixed2) {
            aid_a = ctx.deg1_neighbor(aid1)?;
            aid_b = ctx.closest_neighbor(aid_a, aid2);
            typ = 2;
        } else {
            aid_b = ctx.deg1_neighbor(aid2)?;
            aid_a = ctx.closest_neighbor(aid_b, aid1);
            typ = 3;
        }
        let v2 = self.loc(aid1)?.sub(self.loc(aid_a)?);
        let v1 = self.loc(aid_b)?.sub(self.loc(aid_a)?);
        let cross = v1.x * v2.y - v1.y * v2.x;
        self.apply_open_angle_rotation(aid1, aid2, aid_a, aid_b, typ, cross);
        Ok(())
    }

    /// `removeCollisionsBondFlip`.
    fn remove_collisions_bond_flip(&mut self, ctx: &Ctx) -> DResult<()> {
        let mut colls = self.find_collisions(ctx, true);
        let mut done_bonds: BTreeMap<usize, u32> = BTreeMap::new();
        let mut iter = 0;
        while iter < MAX_COLL_ITERS && !colls.is_empty() {
            let ncols = colls.len();
            let c = colls[0];
            let rot_bonds = ctx.rotatable_bonds(c.0, c.1);
            let prev_density = self.total_density();
            for ri in rot_bonds {
                let cur = done_bonds.get(&ri).copied();
                if cur.is_none() || cur.unwrap() < NUM_BONDS_FLIPS {
                    match cur {
                        None => {
                            done_bonds.insert(ri, 1);
                        }
                        Some(v) => {
                            done_bonds.insert(ri, v + 1);
                        }
                    }
                    self.flip_about_bond(ctx, ri, true)?;
                    colls = self.find_collisions(ctx, true);
                    let mut new_density = self.total_density();
                    if colls.len() < ncols {
                        done_bonds.insert(ri, NUM_BONDS_FLIPS);
                        break;
                    } else if colls.len() == ncols && new_density < prev_density {
                        break;
                    } else {
                        self.flip_about_bond(ctx, ri, true)?;
                        // Recomputed for its density side effect, as RDKit does.
                        self.find_collisions(ctx, true);
                        self.flip_about_bond(ctx, ri, false)?;
                        colls = self.find_collisions(ctx, true);
                        new_density = self.total_density();
                        if colls.len() < ncols {
                            done_bonds.insert(ri, NUM_BONDS_FLIPS);
                            break;
                        } else if colls.len() == ncols && new_density < prev_density {
                            break;
                        } else {
                            self.flip_about_bond(ctx, ri, false)?;
                            colls = self.find_collisions(ctx, true);
                        }
                    }
                }
            }
            iter += 1;
        }
        Ok(())
    }

    /// `removeCollisionsOpenAngles`.
    fn remove_collisions_open_angles(&mut self, ctx: &Ctx) -> DResult<()> {
        for (a, b) in self.find_collisions(ctx, false) {
            self.open_angles(ctx, a, b)?;
        }
        Ok(())
    }

    fn shorten_terminal_bond(
        &mut self,
        ctx: &Ctx,
        atom: usize,
        degree: usize,
        fixed: bool,
        respect_fixed: bool,
    ) -> DResult<()> {
        if degree != 1 || (respect_fixed && fixed) {
            return Ok(());
        }
        let mut location = self.loc(atom)?;
        let neighbour = ctx.deg1_neighbor(atom)?;
        location = location.sub(self.at(neighbour).loc);
        location = location.scale(0.9);
        if location.length() > 0.75 {
            location = location.add(self.at(neighbour).loc);
            self.at(atom).loc = location;
        }
        Ok(())
    }

    /// `removeCollisionsShortenBonds`.
    fn remove_collisions_shorten_bonds(&mut self, ctx: &Ctx) -> DResult<()> {
        let mut colls = self.find_collisions(ctx, false);
        let mut ncols = colls.len();
        let mut iter = 0;
        while ncols > 0 && iter < MAX_COLL_ITERS {
            let (mut aid1, mut aid2) = colls[0];
            let mut fixed1 = self.get(aid1)?.fixed;
            let mut fixed2 = self.get(aid2)?.fixed;
            if fixed1 && fixed2 {
                colls.remove(0);
                ncols = colls.len();
                iter += 1;
                continue;
            }
            let mut deg1 = ctx.degree(aid1);
            let mut deg2 = ctx.degree(aid2);
            if fixed1 || (deg2 > deg1 && !fixed2) {
                std::mem::swap(&mut deg1, &mut deg2);
                std::mem::swap(&mut aid1, &mut aid2);
                std::mem::swap(&mut fixed1, &mut fixed2);
            }
            let _ = fixed1;
            let mut path = ctx.shortest_path(aid1, aid2);
            if path.is_empty() {
                colls.remove(0);
            } else {
                path.remove(0);
                let n_open = ctx.any_non_ring_bonds(aid1, &path)?;
                if n_open > 0 {
                    self.shorten_terminal_bond(ctx, aid1, deg1, fixed1, false)?;
                    self.shorten_terminal_bond(ctx, aid2, deg2, fixed2, true)?;
                } else {
                    let mut rpath = Vec::new();
                    let mut nbr_map = BTreeMap::new();
                    ctx.recurse_deg_two_ring_atoms(aid1, &mut rpath, &mut nbr_map);
                    if rpath.is_empty() {
                        ctx.recurse_deg_two_ring_atoms(aid2, &mut rpath, &mut nbr_map);
                    }
                    let mut move_map: BTreeMap<usize, P> = BTreeMap::new();
                    for &rpi in &rpath {
                        if self.get(rpi)?.fixed {
                            continue;
                        }
                        let nb = nbr_map[&rpi].clone();
                        let mut mv = self.at(nb[0]).loc;
                        mv = mv.add(self.at(nb[1]).loc);
                        mv = mv.scale(0.5);
                        mv = mv.sub(self.loc(rpi)?);
                        mv.normalize()?;
                        mv = mv.scale(COLLISION_THRES);
                        move_map.insert(rpi, mv);
                    }
                    for &rpi in &rpath {
                        let mv = *move_map.entry(rpi).or_default();
                        let a = self.at(rpi);
                        a.loc = a.loc.add(mv);
                    }
                }
                colls = self.find_collisions(ctx, false);
            }
            ncols = colls.len();
            iter += 1;
        }
        Ok(())
    }
}
