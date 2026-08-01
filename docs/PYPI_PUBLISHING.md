# PyPI publishing for ChiralFold

**Status:** `chiralfold==3.6.0` is live on PyPI:
https://pypi.org/project/chiralfold/

```bash
pip install chiralfold==3.6.0
```

## Future releases

Use Trusted Publishing (OIDC) via `.github/workflows/publish-pypi.yml`.

Configure at https://pypi.org/manage/project/chiralfold/settings/publishing/:

| Field | Value |
|-------|-------|
| PyPI project name | `chiralfold` |
| Owner | `Tommaso-R-Marena` |
| Repository name | `ChiralFold` |
| Workflow name | `publish-pypi.yml` |
| Environment name | `pypi` |

**Fork safety:** `publish-pypi.yml` and the `publish-pypi` job in `release.yml` run only when `github.repository == 'Tommaso-R-Marena/ChiralFold'`. Reviewer forks skip publish cleanly (no OIDC failure).

Fallback: set GitHub Actions secret `PYPI_API_TOKEN` to a project-scoped `pypi-...` token.

If you see `invalid-publisher`, the publisher fields do not match the workflow claims — see the troubleshooting table historically used for the first upload of v3.6.0.
