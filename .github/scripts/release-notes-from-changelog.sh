#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 1 || $# -gt 2 ]]; then
  echo "usage: $0 <version-or-tag> [changelog-file]" >&2
  exit 2
fi

version="${1#v}"
changelog_file="${2:-CHANGELOG.md}"
release_heading="## [${version}] - "

if [[ ! -f "$changelog_file" ]]; then
  echo "changelog not found: $changelog_file" >&2
  exit 1
fi

awk \
  -v release_heading="$release_heading" \
  -v version="$version" \
  -v changelog_file="$changelog_file" '
  !capturing && index($0, release_heading) == 1 {
    capturing = 1
  }
  capturing && index($0, "## [") == 1 && index($0, release_heading) != 1 {
    exit
  }
  capturing {
    print
    if ($0 !~ /^[[:space:]]*$/) {
      nonblank_lines++
    }
  }
  END {
    if (!capturing || nonblank_lines < 2) {
      printf "no release notes found for v%s in %s\n", version, changelog_file > "/dev/stderr"
      exit 1
    }
  }
' "$changelog_file"
