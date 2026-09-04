# Notes for outputs

All runs were performed using the parameters in example `config.json` and/or specified in relevant .md files.

Targets are considered "bad", "contaminated" and/or "broken" and are excluded if they (for all ObsIDs):
- Do not contain the minimum 30 counts.
- Fail to complete a run in the pipeline.
- Are contaminated by other sources.
- Errors are pinned at boundaries (0.0, 5.0).
- Net counts at end of pipeline ≤ 0.0.
  
| Type | Total | Test | Normal | Weak | Strong |
|---|---|---|---|---|---|
| Chosen | 91 | 6 | 44 | 31 | 10 |
| Successful* | 32 | 6 | 19 | 6 | 1 |
| Percentage | 35.1% | 100% | 43.2% | 19.4% | 10% |

<sub><sup>* Successful runs are only ones that succeed the parameters described. They do not exclude those manually removed or are flagged in `examples.md. </sub></sup>

Specific folder shows outputs produced by each target for that target after complete run. Only showing specific ones for ease of use.
General folder shows outputs produced regardless of target, some depend on target_catalogue and previous runs/existing data.
