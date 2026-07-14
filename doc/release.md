# GitHub release procedure

This project uses Cargo release-mode builds to produce binaries for GitHub
Releases. It is not published to crates.io. Pushing the release tag triggers
the workflow in `.github/workflows/release.yml`.

## Prepare the release candidate

1. Start from the intended release commit with no uncommitted changes and
   confirm the normal CI workflow is green.
2. Set the package version in `Cargo.toml` and `Cargo.lock`, then give the
   matching changelog section its release date.
3. Run the local release checks:

   ```bash
   rustup toolchain install 1.97.0 --profile minimal \
     --component rustfmt --component clippy
   rustup toolchain install 1.88.0 --profile minimal
   cargo +1.97.0 install cargo-audit --locked --version 0.22.2

   cargo +1.97.0 fmt --all -- --check
   cargo +1.97.0 clippy --locked --all-targets --all-features -- -D warnings
   cargo +1.97.0 test --locked --all-targets --all-features
   RUSTDOCFLAGS="-D warnings -D missing_docs" \
     cargo +1.97.0 doc --locked --no-deps --all-features
   cargo +1.97.0 test --locked --doc --all-features
   cargo +1.88.0 test --locked --lib
   cargo +1.97.0 audit --deny warnings
   cargo +1.97.0 package --locked
   cargo +1.97.0 build --release --locked --bin trackclustertu
   ```

   These versions match the tagged-release workflow. Installing an unpinned
   `cargo-audit` or running the unqualified default toolchain does not reproduce
   the release gate.

   `cargo package` is a source-completeness check only. Its `.crate` output is
   not uploaded or published.

4. Verify the local binary and intended tag:

   ```bash
   version="$(cargo +1.97.0 pkgid --locked | sed 's/.*@//')"
   tag="v${version}"
   test "$(./target/release/trackclustertu --version)" = \
     "trackclustertu ${version}"
   ./target/release/trackclustertu --help >/dev/null
   test -z "$(git tag --list "$tag")"
   ```

## Publish the GitHub release

Commit and push the reviewed release candidate before creating the tag. Then
create an annotated tag on that exact commit:

```bash
version="$(cargo +1.97.0 pkgid --locked | sed 's/.*@//')"
tag="v${version}"
git tag -a "$tag" -m "trackclusterTU ${tag}"
test "$(git rev-parse "${tag}^{}")" = "$(git rev-parse HEAD)"
git push origin "$tag"
```

The workflow rejects a tag that does not match the Cargo package version. It
runs the release gate, builds all supported targets, creates one `SHA256SUMS`
manifest, and publishes the GitHub release automatically. Each build receives
the tag commit through `TRACKCLUSTERTU_GIT_REVISION`; the smoke test rejects a
binary that does not contain those exact revision bytes. `run_manifest.json`
reads the same embedded compile-time value, but the release smoke test does not
generate or inspect a manifest.

## Verify the published release

1. Confirm the GitHub release contains these three archives and one
   `SHA256SUMS` manifest:

   - `trackclustertu-<tag>-x86_64-unknown-linux-musl.tar.gz`
   - `trackclustertu-<tag>-aarch64-unknown-linux-gnu.tar.gz`
   - `trackclustertu-<tag>-aarch64-apple-darwin.tar.gz`

2. Download all four assets into an empty directory and verify the archives:

   ```bash
   sha256sum --check SHA256SUMS
   ```

   On macOS, use `shasum -a 256 --check SHA256SUMS` instead.

3. Confirm each archive contains exactly `trackclustertu`, `LICENSE`, and
   `README.md`, then smoke-test each binary on its native platform:

   ```bash
   tar -tzf trackclustertu-<tag>-<target>.tar.gz | sort
   test "$(./trackclustertu --version)" = "trackclustertu <version>"
   ./trackclustertu --help >/dev/null
   ```

4. Record the release tag commit, `SHA256SUMS`, workflow result, and native
   smoke-test results in the release record maintained outside this repository.
