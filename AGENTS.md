# AGENTS.md — How to work in this repository

## TL;DR
- You are contributing to **AGAT** (Perl).
- Default to incremental changes + tests. Prefer small PRs.
- It is a suite of multiple scripts. Some contain similar problems, check all scripts before fixing in only one of them.

## Agent quickstart (isolated shells; no sudo)

```bash
- Bootstrap once: `bash .agents/bootstrap.sh`
- Activate env per step: use `.agents/with-perl-local.sh`
- (smoke tests are not ready now) Run **smoke** only while iterating: `prove -lr t/smoke`
- Run **full test suite** only before a commit/PR: `make test`
```

## Environment
- Perl versions: 5.36–5.42 preferred.

## Tests
- Full suite: `make test`
- Coverage (optional): `cover -test` (uploader: Coveralls/Codecov if configured)
- Config management: `agat config --expose 2>&1 1>/dev/null` (redirecting to suppress messages). By default, it creates `agat_config.yaml` in the workdir. Use `--output <local_config_path>` to rename and `--config <local_config_path>` when invoking scripts.

### Default parser config
If not overriden by --expose invokation with custom parameters, tests use these `share/agat_config.yaml` defaults:
```yaml
output_format: GFF            # output format
gff_output_version: 3         # GFF version
gtf_output_version: relax     # accept all feature types
deflate_attribute: false      # keep multi-value attrs
verbose: 1                    # verbosity 0–4
progress_bar: true            # show progress bar
log: true                     # write log file
debug: false                  # extra verbosity
tabix: false                  # no tabix output
merge_loci: false             # don't merge overlapping loci
throw_fasta: false            # keep embedded FASTA
force_gff_input_version: 0    # auto-detect input version
create_l3_for_l2_orphan: true # add exon if missing
locus_tag: [locus_tag, gene_id] # fallback attributes
prefix_new_id: agat           # ID prefix
clean_attributes_from_template: false # keep attrs
check_sequential: true
check_l2_linked_to_l3: true
check_l1_linked_to_l2: true
remove_orphan_l1: true
check_all_level3_locations: true
check_cds: true
check_exons: true
check_utrs: true
check_all_level2_locations: true
check_all_level1_locations: true
check_identical_isoforms: true
```

### Test data
- Use fixtures in `t/data/`. Do not download external corpora.
- Max sample size ~ a few MB; keep runtime deterministic (no network).

## Style & quality
- Follow existing patterns (Moose, strict/warnings). New modules go in `lib/AGAT/...`.
- Prefer pure-Perl core and explicit modules; document non-core deps in `Makefile.PL` and `cpanfile` (if present).
- Lint (if `perlcritic` is installed): `perlcritic --gentle lib bin`
  - Install via `cpm install Perl::Critic` (or `cpanm Perl::Critic`), otherwise skip

## Docs
- User docs live under `docs/` (MkDocs). Preview with: `mkdocs build`
  - Requires `mkdocs`; install with `pip install mkdocs` or skip if unavailable
- If a tool’s CLI changes, update its `.md` page and usage examples.

## CI mirror
- CI runs the same build/test commands on Linux with `actions-setup-perl`.
- Keep PRs green by running the **fast suite** locally before pushing.

## Boundaries
- Don’t edit vendored/generated files or CI YAML except as requested.
- Avoid O(n²) over whole-genome GFFs; stream where practical.

## PR checklist
- [ ] Unit tests for new behavior
- [ ] Updated docs/examples
- [ ] Changelog entry (if required)
- [ ] CI green on Linux; Windows/macOS issues noted if relevant
