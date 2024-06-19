# Run the R scripts

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

Run Script directly in R:
```
system("Rscript /Users/ju5263ta/Github/Monocytes/Scripts/240607_Controls_ACTb.R")
system("Rscript /Users/ju5263ta/Github/Monocytes/Scripts/240607_Controls_B2M.R")
system("Rscript /Users/ju5263ta/Github/Monocytes/Scripts/240607_Stroke_B2M.R")
```

All output will be saved in a folder: TODAYS-DATE_Stroke/Controls_HK_Results