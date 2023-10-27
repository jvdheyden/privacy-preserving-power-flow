# Instructions for replicating our results

We assume basic familiarity with MP-SPDZ. First, move the .mpc files into
`MP-SPDZ/Programs/Source/` and move the data folder as well as `lib.py` into
the MP-SPDZ top level directory. Create player data and setup SSL as required by
MP-SPDZ. To compile the programs for the rural example, make sure that within 
the .mpc files the import is set to data.rural. For the semi-urban example,
the import should be set to data.simbench. To accommodate the precision specified
in our .mpc files, add the following to `MP-SPDZ/CONFIG.mine` and remake.
```
MY_CFLAGS += -DINSECURE                                                         
MOD = -DRING_SIZE=261 -DGFP_MOD_SZ=3 -DMAX_MOD_SZ=12
```

You can now compile the .mpc files using the following commands:
```
./compile.py -Y -I -O newton_raphson_lu; ./compile.py -Y -I -O newton_raphson_GMRES
``` 

To generate correlated randomness for offline preprocessing for 13 players (rural grid), use this command:
```
./Fake-Offline.x 13 -lgp 192 -e 1,32,33,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,62,63,64,65,66,67,72,97,102,104,105,107,126,127,128,129,130,131,138,166,168,169,171,190,191,192,194,227,230,232,234,235,265,362
```

For 39 players (semi-urban grid), use this command:

```
./Fake-Offline.x 39 -lgp 192 -e 1,25,26,27,28,29,30,31,32,33,36,37,38,39,40,41,42,43,44,45,46,47,48,49,52,62,63,64,65,66,67,72,98,102,104,105,107,126,127,128,129,130,131,138,144,166,168,169,171,194,234,320
```

Now run benchmarks for various protocols, e.g. Shamir:

```
PLAYERS=13 Scripts/shamir.sh -v -F newton_raphson_lu
```

To run offline and online phase together, remove the -F flag. To compile the .mpc files for ring-
based protocols, use the following commands:
```
./compile.py -Y -I -O -R 261 newton_raphson_lu; ./compile.py -Y -I -O -R 261 newton_raphson_GMRES
```
And to generate correlated randomness for ring-based protocols:
```
./Fake-Offline.x 13 -Z 261 -e 1,36,37,38,39,40,41,42,43,44,45,46,47,48,49,52,64,65,104,126,127,128,129,130,131,166,168,169,171,194,234
```
resp.

```
./Fake-Offline.x 39 -Z 261 -e 1,25,26,27,28,29,30,31,32,33,36,37,38,39,40,41,42,43,44,45,46,47,48,49,52,62,63,64,65,66,67,72,98,102,104,105,107,126,127,128,129,130,131,138,144,166,168,169,171,194,234,320
```
To simulate a certain RTT between prosumers you can use Linux's `tc`:
```
sudo tc  qdisc add dev lo root netem delay 10ms
```

All benchmarks were produced using the MP-SPZ version at commit b35d0d5.
