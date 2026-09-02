# Release-key custody

Release artifacts must be signed by a key held outside the repository. The
private key must never be committed, uploaded as an artifact, printed in logs,
or placed in a normal repository variable. Store it as the masked GitHub
Actions secret `SCHEMATIC_RELEASE_PRIVATE_KEY` (or in an equivalent hardware
backed release system) and restrict its use to the protected release
environment.

The matching public key is published separately with the release evidence and
its SHA-256 fingerprint is recorded in the release notes. Consumers verify the
exact provenance bytes with:

```sh
scripts/verify_release_key.sh \
  schematic.provenance.json schematic.provenance.sig \
  schematic-release-public-key.pem <published-fingerprint>
```

The verifier checks both the DER-encoded public-key fingerprint and the
detached SHA-256 signature. A key rotation requires a new fingerprint, an
explicit release note, and verification of the first artifact with both the
old and new operational records. A local disposable key is test evidence only
and must not be presented as release provenance.

## Activation checklist

- [ ] Generate the production key in the maintainer-controlled secret store.
- [ ] Register `SCHEMATIC_RELEASE_PRIVATE_KEY` in the protected GitHub
  environment; do not expose it to pull requests.
- [ ] Publish the matching public key and fingerprint through the project
  release evidence channel.
- [ ] Run the manual `Release key evidence` workflow and retain its artifact.
- [ ] Verify the signed provenance from a clean checkout using the published
  public key.
