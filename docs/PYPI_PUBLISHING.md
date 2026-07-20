# PyPI publishing for ChiralFold

`pip install chiralfold` requires a successful upload to
https://pypi.org/project/chiralfold/.

## Current failure mode

GitHub Actions logs show:

```text
Trusted publishing exchange failure:
* invalid-publisher: valid token, but no corresponding publisher
  (Publisher with matching claims was not found)
```

OIDC works on the GitHub side; **PyPI has no Trusted Publisher that matches
this repository's claims**. That is a PyPI project setting, not a code bug.

## Fix option A — Trusted Publisher (preferred)

1. Log in to PyPI as the owner of project `chiralfold`
   (or create the project on first publish).
2. Open **Project settings → Publishing → Add a new pending publisher**.
3. Enter **exactly**:

| Field | Value |
|-------|-------|
| PyPI project name | `chiralfold` |
| Owner | `Tommaso-R-Marena` |
| Repository name | `ChiralFold` |
| Workflow name | `publish-pypi.yml` |
| Environment name | `pypi` |

4. Ensure the GitHub repo has an Environment named **`pypi`**
   (Settings → Environments). No secrets are required for OIDC.
5. Re-run **Publish to PyPI** (Actions → workflow_dispatch) or publish a
   GitHub Release / tag `v*`.

Claims emitted by the failing run (for comparison):

- `repository`: `Tommaso-R-Marena/ChiralFold`
- `workflow_ref`: `.../.github/workflows/publish-pypi.yml@...`
- `environment`: `pypi`

If you configure the publisher **without** an Environment name, either leave
the Environment field blank on PyPI **and** remove `environment: pypi` from
`.github/workflows/publish-pypi.yml`, or keep both sides set to `pypi`.

## Fix option B — API token fallback

1. Create a PyPI API token scoped to project `chiralfold`
   (Account settings → API tokens).
2. Add GitHub Actions secret **`PYPI_API_TOKEN`** with that token value
   (`pypi-...`).
3. Re-run the publish workflow. The workflow passes
   `password: ${{ secrets.PYPI_API_TOKEN }}` into
   `pypa/gh-action-pypi-publish`; when the secret is set, token auth is used
   instead of OIDC.

## Verify

```bash
pip index versions chiralfold
# or
pip install chiralfold==3.5.1
python -c "import chiralfold; print(chiralfold.__version__)"
```

## Related workflows

- `.github/workflows/publish-pypi.yml` — manual / release publish
- `.github/workflows/release.yml` — tag `v*` release (also publishes when
  Trusted Publisher or `PYPI_API_TOKEN` is configured)
