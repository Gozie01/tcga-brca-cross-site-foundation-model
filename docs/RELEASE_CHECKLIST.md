# Release checklist

- [ ] Confirm author names, order, affiliations, and ORCIDs.
- [ ] Confirm repository title and manuscript title.
- [ ] Scan for credentials, access tokens, email addresses not intended for publication, and local absolute paths.
- [ ] Confirm no raw WSI, controlled-access data, UNI weights, or third-party papers are included.
- [ ] Validate all CSV and JSON files and the NumPy feature matrix.
- [ ] Run the analysis smoke tests and figure synthesis.
- [ ] Regenerate `checksums_sha256.txt` after the final change.
- [ ] Commit and push release contents to GitHub.
- [ ] Enable the repository in Zenodo.
- [ ] Create GitHub release `v1.0.0`.
- [ ] Verify the Zenodo software record and DOI.
- [ ] Create and cross-link a separate Zenodo dataset record if the derived outputs are deposited separately.
- [ ] Insert final DOI(s) into the manuscript Data and Code Availability statements.
