# pyGecko

The main goal of this project is to build an modular library for GC-MS and GC-FID data analysis that can be easily extended and maintained.

## 1. Think Before Coding

**Don't assume. Don't hide confusion. Surface tradeoffs.**

Before implementing:
- State your assumptions explicitly. If uncertain, ask.
- If multiple interpretations exist, present them - don't pick silently.
- If a simpler approach exists, say so. Push back when warranted.
- If something is unclear, stop. Name what's confusing. Ask.

## 2. Simplicity First

**Minimum code that solves the problem. Nothing speculative.**

- No features beyond what was asked.
- No abstractions for single-use code.
- No "flexibility" or "configurability" that wasn't requested.
- No error handling for impossible scenarios.
- If you write 200 lines and it could be 50, rewrite it.

Ask yourself: "Would a senior engineer say this is overcomplicated?" If yes, simplify.

## 3. Surgical Changes

**Touch only what you must. Clean up only your own mess.**

When editing existing code:
- Don't "improve" adjacent code, comments, or formatting.
- Don't refactor things that aren't broken.
- Match existing style, even if you'd do it differently.
- If you notice unrelated dead code, mention it - don't delete it.

When your changes create orphans:
- Remove imports/variables/functions that YOUR changes made unused.
- Don't remove pre-existing dead code unless asked.

The test: Every changed line should trace directly to the user's request.

## 4. Goal-Driven Execution

**Define success criteria. Loop until verified.**

Transform tasks into verifiable goals:
- "Add validation" → "Write tests for invalid inputs, then make them pass"
- "Fix the bug" → "Write a test that reproduces it, then make it pass"
- "Refactor X" → "Ensure tests pass before and after"

For multi-step tasks, state a brief plan:
```
1. [Step] → verify: [check]
2. [Step] → verify: [check]
3. [Step] → verify: [check]
```

Strong success criteria let you loop independently. Weak criteria ("make it work") require constant clarification.


## 5. Testing

Use the `python-testing-patterns` skill for all test work. Consult it before writing tests, designing a test suite, or changing testing infrastructure.

Follow the TDD cycle without exception:

1. **RED** — write a failing test that specifies the desired behavior. Run it and confirm it fails for the expected reason.
2. **GREEN** — write the minimal implementation that makes the test pass. No extra features, no speculative abstraction.
3. **REFACTOR** — clean up implementation and test code while keeping the suite green.

Do not write implementation code before a failing test exists for it.

### Requirements

- pytest only. Use `pytest.raises` for exceptions, fixtures for setup/teardown, `@pytest.mark.parametrize` for input variation.
- Minimum 80% coverage overall; 100% on critical paths. Verify with `pytest --cov=pyGecko --cov-report=term-missing`.
- Mock external dependencies (network, database, filesystem). Tests must run offline and be independent of execution order.
- Mark slow and integration tests so `pytest -m "not slow"` stays fast.
- Test behavior, not implementation details.

## 6. Documentation Roles

Three documents serve distinct purposes — keep them separate and use each for its intended role:

- **AGENTS.md** (this file; `CLAUDE.md` is a symlink to it, so editing either one edits both) — dictates
  **general agent behaviour**: how to think, code, and work in this repo. It is about *conduct*, not
  architecture or usage.
- **docs/architecture.md** — used to **plan the architecture to implement and document design decisions**.
  Update it *after* changes have been implemented so it reflects the current design and the reasoning
  behind it.
- **README.md** — **outward-facing**: describes how to use the agent and how to contribute to it.
