---
name: always-push-commits
description: "User wants every commit pushed to the remote immediately, without being asked"
metadata: 
  node_type: memory
  type: feedback
  originSessionId: 8fbc0ac7-489a-4716-ba8a-72a7aa11fee1
---

When I make a commit, push it to the remote (`git push`) without waiting to be asked.
At the **end of every work batch**, commit and push the work — don't leave a finished
batch uncommitted in the working tree.

**Why:** The user stated "always push commits" (2026-06-10) and reaffirmed "always
commit and push work at the end of a work batch" (2026-06-10). They treat a commit as
not-done until it's on the remote, and a batch of work as not-done until it's committed
*and* pushed; leaving work local is an extra round-trip they don't want.

**How to apply:** When a coherent batch of work is finished, `git commit` it and then
`git push` to the tracking branch — without being asked. This overrides the default
"commit/push only when asked" stance for this project. Still branch/PR per the repo's
convention where one exists — here history is linear on `main`, so push straight to
`main`. See [[flatten-repo-root-gotcha]].
