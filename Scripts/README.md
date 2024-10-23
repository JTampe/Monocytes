# Run the R scripts as follows

set up the Terminal:
```
cd Github/Monocytes
git remote set-url origin git@github.com:JTampe/Monocytes.git

```

Start R in the Data folder:
```
cd Data 
R
system( "ls" )
```

Run Script for the Control Data directly in R via the terminal:
```
system("Rscript /Users/ju5263ta/Github/Monocytes/Scripts/240905_Controls_MFIregression.R")
```

Run Script for the Stroke Data directly in R via the terminal:
```
system("Rscript /Users/ju5263ta/Github/Monocytes/Scripts/240920_Stroke_NormAll_imp.R")
```



All output will be saved in a folder: TODAYS-DATE_Stroke/Controls_Github_Results
