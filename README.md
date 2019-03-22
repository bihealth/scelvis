# single-cell visualization in SODAR using dash + docker

prototype still based on locally provided datasets 

## conda

the dash environment is described in `dash_env.hml` and should install everything that's needed

## first trial version

run `python hnc/hnc.py` and then go to http://127.0.0.1:8050 to have a look

## first version of full app

the app reads the variable `DASH_DATADIR` and then uses data from that directory

```
cd app
export DASH_DATADIR=datasets/hgmm_1k
python app.py
```

