# Results of Limited CCD0 Expansion

## EBOV (516)

| Method               | log(MCC)    | Runtime |
| -------------------- | ----------- | ------- |
| CCD1                 | -1752.89    | 171s    |
| HIPSTR               | −950.08     | **57s** |
| CCD0: no expansion   | -950.08     | 182s    |
| CCD0: 50n expansion  | -946.77     | 171s    |
| CCD0: 100n expansion | -939.77     | 176s    |
| CCD0                 | **-923.66** | 17 min  |

## EBOV (1610)

| Method               | log(MCC)     | Runtime |
| -------------------- | ------------ | ------- |
| CCD1                 | -4686.60     | 22s     |
| HIPSTR               | −3039.60     | **6s**  |
| CCD0: no expansion   | -3039.60     | 25s     |
| CCD0: 50n expansion  | -2931.15     | 22s     |
| CCD0: 100n expansion | -2921.48     | 26s     |
| CCD0                 | **-2886.24** | 64s     |

## Simulated EBOV (1610)

| Method               | log(MCC)     | Runtime |
| -------------------- | ------------ | ------- |
| CCD1                 | -3543.36     | 194s    |
| HIPSTR               | -1786.64     | **81s** |
| CCD0: no expansion   | -1786.64     | 162s    |
| CCD0: 50n expansion  | -1775.41     | 184s    |
| CCD0: 100n expansion | -1781.02     | 204s    |
| CCD0                 | **-1767.57** | 17min   |

## SARS-CoV-2 (3959)

| Method               | log(MCC)     | Runtime |
| -------------------- | ------------ | ------- |
| CCD1                 | -11540.45    | 65s     |
| HIPSTR               | -8630.85     | **21s** |
| CCD0: no expansion   | -8745.53     | 57s     |
| CCD0: 50n expansion  | -8470.35     | 86s     |
| CCD0: 100n expansion | -8448.67     | 201s    |
| CCD0                 | **-8416.94** | 11min   |

## SARS-CoV-2 (15616)

| Method               | log(MCC)      | Runtime |
| -------------------- | ------------- | ------- |
| CCD1                 | -58636.38     | 150s    |
| HIPSTR               | -57546.68     | **85s** |
| CCD0: no expansion   | -51399.93     | 144s    |
| CCD0: 50n expansion  | -50785.12     | 18min   |
| CCD0: 100n expansion | -50458.26     | 77min   |
| CCD0                 | **-49922.53** | 9h      |
