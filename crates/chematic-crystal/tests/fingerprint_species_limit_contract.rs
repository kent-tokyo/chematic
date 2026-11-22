use chematic_core::Element;
use chematic_crystal::fingerprint::{FingerprintError, FingerprintLimits, geometry_fingerprint};
use chematic_crystal::{
    FractionalCoord, Lattice, Occupancy, PeriodicSite, PeriodicStructure, SiteSpecies,
};

#[test]
fn mixed_site_species_limit_counts_species_independently_of_sites() {
    let site = PeriodicSite::new(
        vec![
            SiteSpecies {
                element: Element::FE,
                occupancy: Occupancy::new(0.5).unwrap(),
            },
            SiteSpecies {
                element: Element::NI,
                occupancy: Occupancy::new(0.5).unwrap(),
            },
        ],
        FractionalCoord::new([0.0, 0.0, 0.0]),
        None,
    )
    .unwrap();
    let structure = PeriodicStructure::new(Lattice::cubic(4.0).unwrap(), vec![site]).unwrap();
    let limits = FingerprintLimits {
        max_sites: 1,
        max_species: 2,
    };
    let result = geometry_fingerprint(Some(&structure), "mixed", limits).unwrap();
    assert_eq!((result.site_count, result.species_count), (1, 2));
    let restricted = FingerprintLimits {
        max_species: 1,
        ..limits
    };
    assert_eq!(
        geometry_fingerprint(Some(&structure), "mixed", restricted),
        Err(FingerprintError::ResourceLimit)
    );
}
