# PyPI publishing for ChiralFold

**Status:** `chiralfold==3.5.1` is live on PyPI:
https://pypi.org/project/chiralfold/

```bash
pip install chiralfold==3.5.1
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

Fallback: set GitHub Actions secret `PYPI_API_TOKEN` to a project-scoped `pypi-...` token.

If you see `invalid-publisher`, the publisher fields do not match the workflow claims — see the troubleshooting table historically used for the first upload of v3.5.1.
