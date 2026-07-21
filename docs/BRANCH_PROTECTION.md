# Branch protection for master

Jobs in `.github/workflows/ci.yml` already run **in parallel**. A final
**`CI success`** job waits for every sibling job and fails if any did.

## Required status check (blocks merge until green)

In GitHub:

1. **Settings → Branches → Add/Edit branch protection rule** for `master`
2. Enable **Require a pull request before merging**
3. Enable **Require status checks to pass before merging**
4. Enable **Require branches to be up to date before merging** (recommended)
5. Search and select the single check: **`CI success`**
6. Optionally also require **`HF Space build + smoke`** from
   `deploy-hf-space.yml` if you want Space Docker smoke on every PR that
   touches `hf_space/` / `web/` / `chiralfold/`

Do **not** list every matrix cell (Linux 3.9, 3.10, …) as separate required
checks — the `CI success` gate already aggregates them.

## Why one gate

- Matrix jobs rename / get added without updating branch rules
- A single green check is clearer for reviewers
- `if: always()` on the gate still runs when an upstream job fails, so the
  PR shows a definitive fail rather than a skipped required check
