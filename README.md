# Prerequisites

* Create file `.signals/auth.json` with secrets from
  `lu/signals-sandbox API key and LuPySNB setup` item from Keeper

```json
{
  "base_url": "...",
  "api_key": ""
}
```

# Setup

1. `poetry install` to install dependencies
2. `poetry run python forward.py` to run "offline" forward flow on
   examples from SNB:Exp-134.  Output will be saved in `data/forward`
   1. manually uncomment the desired line `plate_container_eid =
      'plateContainer:<eid>'` in `if __name__ == '__main__'` to run a
      specific plate container from experiment
3. `poetry run python return.py` to run "offline" return flow pipeline
   that converts output XLSX to IDS-like JSON
