# Private-to-public publication checklist

This starter repository is a **local deliverable**. No remote GitHub repository has been created and no license has been chosen.

1. Confirm that laboratory/institutional policies and coauthors permit sharing these Python scripts. Get approval of any licensing terms and attribution.
2. Review `git status` and `git diff --cached` before committing. Do not commit measured raw spectra, confidential metadata, unpublished results, personal file paths, passwords, or other non-shareable files.
3. Run `python -m pytest -q` in a fresh environment and compare at least one representative approved measurement against the known original script outputs.
4. Test the README installation and command examples on at least one Windows lab computer.
5. Optionally publish a **private GitHub repository** first, invite collaborators, and request scientific review.
6. Only make the repository public after ownership, code/data sharing, and licensing review.

## Suggested Git commands (after review)

Run from the unzipped repository root after installing Git:

```bash
git init
git add .
git status
git commit -m "Initial Raman processing tools with examples and tests"
```

Create an empty **private** repository under the approved GitHub account/organization, then follow GitHub's instructions to add its remote and push. Avoid hardcoding a personal account or organization in these instructions.
