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
   cargo fmt --all -- --check
   cargo clippy --locked --all-targets --all-features -- -D warnings
   cargo test --locked --all-targets --all-features
   RUSTDOCFLAGS="-D warnings -D missing_docs" \
     cargo doc --locked --no-deps --all-features
   cargo test --locked --doc --all-features
   cargo +1.88.0 test --locked --lib
   cargo audit --deny warnings
   cargo package --locked
   cargo build --release --locked --bin trackclustertu
   ```

   `cargo package` is a source-completeness check only. Its `.crate` output is
   not uploaded or published.

4. Verify the local binary and intended tag:

   ```bash
   version="$(cargo pkgid --locked | sed 's/.*@//')"
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
version="$(cargo pkgid --locked | sed 's/.*@//')"
tag="v${version}"
git tag -a "$tag" -m "trackclusterTU ${tag}"
test "$(git rev-parse "${tag}^{}")" = "$(git rev-parse HEAD)"
git push origin "$tag"
```

The workflow rejects a tag that does not match the Cargo package version. It
runs the release gate, builds all supported targets, creates one `SHA256SUMS`
manifest, and publishes the GitHub release automatically. Each build receives
the tag commit through `TRACKCLUSTERTU_GIT_REVISION`; the smoke test rejects a
binary that does not embed that exact revision for `run_manifest.json`.

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
