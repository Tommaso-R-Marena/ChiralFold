# PyPI publishing for ChiralFold

`pip install chiralfold` requires a successful upload to
https://pypi.org/project/chiralfold/.

## Diagnosis from the latest failing run

Build of `chiralfold-3.5.1` **succeeds**. Upload fails with:

```text
Trusted publishing exchange failure:
* invalid-publisher: valid token, but no corresponding publisher
  (Publisher with matching claims was not found)
```

GitHub OIDC is working. **PyPI has no Trusted Publisher (or pending publisher)
that matches these claims** from the failing job:

| Claim | Value from the log |
|-------|--------------------|
| `repository` | `Tommaso-R-Marena/ChiralFold` |
| `workflow_ref` | `.../.github/workflows/publish-pypi.yml@refs/heads/master` |
| `environment` | `pypi` |
| `sub` | `repo:Tommaso-R-Marena/ChiralFold:environment:pypi` |

This is a PyPI account/project setting, not a packaging bug.

## Fix option A — Pending / Trusted Publisher (preferred)

Because `chiralfold` is not on PyPI yet, use a **pending publisher**:

1. Log in to PyPI as **Tommaso-R-Marena** (or the account that should own the project).
2. Open **Your account → Publishing → Add a new pending publisher**
   (https://pypi.org/manage/account/publishing/).
3. Enter **exactly** (case-sensitive):

| Field | Value |
|-------|-------|
| PyPI project name | `chiralfold` |
| Owner | `Tommaso-R-Marena` |
| Repository name | `ChiralFold` |
| Workflow name | `publish-pypi.yml` |
| Environment name | `pypi` |

4. Confirm GitHub has an Environment named **`pypi`**
   (repo → Settings → Environments). No secrets required for OIDC.
5. Re-run **Actions → Publish to PyPI → Run workflow** (leave tag blank to
   publish current `master`, or set `v3.5.1`).

Common mismatches that cause `invalid-publisher`:

- Workflow name includes a path (`/.github/workflows/publish-pypi.yml`) — use **only** `publish-pypi.yml`
- Environment left blank on PyPI while the workflow has `environment: pypi` (or the reverse)
- Repo owner typo / wrong capitalization (`tommaso-r-marena` vs `Tommaso-R-Marena` — use the GitHub owner string PyPI expects; usually the login as shown on GitHub)
- Publisher registered under a different PyPI account than the one you expect

## Fix option B — API token fallback

1. Create a PyPI API token (Account settings → API tokens). For a new project,
   use a **user-scoped** token for the first upload, then switch to a
   project-scoped token.
2. Add GitHub Actions secret **`PYPI_API_TOKEN`** with the token value (`pypi-...`).
3. Re-run the publish workflow. When that secret is set, token auth is used
   instead of OIDC.

## Verify

```bash
pip index versions chiralfold
pip install chiralfold==3.5.1
python -c "import chiralfold; print(chiralfold.__version__)"
```

## Related workflows

- `.github/workflows/publish-pypi.yml` — manual / release publish
- `.github/workflows/release.yml` — tag `v*` release
