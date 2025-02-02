# Results of Limited CCD0 Expansion

## EBOV (516)

| Method               | log(MCC)    | Runtime |
| -------------------- | ----------- | ------- |
| CCD1                 | -1752.89    | 171s    |
| HIPSTR               | −950.08     | **57s** |
| CCD0: no expansion   | -950.08     | 197s    |
| CCD0: 50n expansion  | -946.77     | 168s    |
| CCD0: 100n expansion | -939.77     | 148s    |
| CCD0                 | **-923.66** | 11 min  |

## EBOV (1610)

| Method               | log(MCC)     | Runtime |
| -------------------- | ------------ | ------- |
| CCD1                 | -4686.60     | 22s     |
| HIPSTR               | −3039.60     | **6s**  |
| CCD0: no expansion   | -3039.60     | 23s     |
| CCD0: 50n expansion  | -2931.15     | 25s     |
| CCD0: 100n expansion | -2896.87     | 32s     |
| CCD0                 | **-2864.92** | 44s     |

## Simulated EBOV (1610)

| Method               | log(MCC)     | Runtime |
| -------------------- | ------------ | ------- |
| CCD1                 | -3543.36     | 194s    |
| HIPSTR               | -1786.64     | **81s** |
| CCD0: no expansion   | -1786.64     | 228s    |
| CCD0: 50n expansion  | 1775.41      | 204s    |
| CCD0: 100n expansion | -1774.31     | 204s    |
| CCD0                 | **-1760.87** | 6min    |

## SARS-CoV-2 (3959)

| Method               | log(MCC)     | Runtime |
| -------------------- | ------------ | ------- |
| CCD1                 | -11540.45    | 65s     |
| HIPSTR               | -8630.85     | **21s** |
| CCD0: no expansion   | -8630.85     | 57s     |
| CCD0: 50n expansion  | -8470.35     | 86s     |
| CCD0: 100n expansion | -8327.63     | 147s    |
| CCD0                 | **-8325.58** | 6min   |

## SARS-CoV-2 (15616)

| Method               | log(MCC)      | Runtime |
| -------------------- | ------------- | ------- |
| CCD1                 | -58636.38     | 150s    |
| HIPSTR               | -51352.61     | **85s** |
| CCD0: no expansion   | -51352.61     | 146s    |
| CCD0: 50n expansion  | -50785.12     | 18min   |
| CCD0: 100n expansion | -50418.37     | 82min   |
| CCD0                 | **-49836.26** | 9h      |
