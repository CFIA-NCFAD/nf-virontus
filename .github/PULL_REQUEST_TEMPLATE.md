<!--
# peterk87/nf-virontus pull request

Many thanks for contributing to peterk87/nf-virontus!

Please fill in the appropriate checklist below (delete whatever is not relevant).
These are the most common things requested on pull requests (PRs).

Remember that PRs should be made against the dev branch, unless you're preparing a pipeline release.

Learn more about contributing: [CONTRIBUTING.md](https://github.com/peterk87/nf-virontus/tree/master/.github/CONTRIBUTING.md)
-->

## PR checklist

- [ ] This comment contains a description of changes (with reason).
- [ ] If you've fixed a bug or added code that should be tested, add tests!
  - [ ] If you've added a new tool - add to the software_versions process and a regex to `scrape_software_versions.py`
  - [ ] If you've added a new tool - have you followed the pipeline conventions in the [contribution docs](https://github.com/peterk87/nf-virontus/tree/master/.github/CONTRIBUTING.md)
  - [ ] If necessary, also make a PR on [the peterk87/nf-test-datasets repo](https://github.com/peterk87/nf-test-datasets/pull/new)
- [ ] Make sure your code lints (`nf-core lint .`).
- [ ] Ensure the test suite passes (`nextflow run . -profile test,docker`).
- [ ] Usage Documentation in `docs/usage.md` is updated.
- [ ] Output Documentation in `docs/output.md` is updated.
- [ ] `CHANGELOG.md` is updated.
- [ ] `README.md` is updated (including new tool citations and authors/contributors).
