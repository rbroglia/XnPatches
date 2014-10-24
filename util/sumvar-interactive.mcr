#!MC 1120

$!VarSet |MFBD| = './'
$!VARSET |ZNUMS|=""
$!VARSET |VNAME|=""
$!VARSET |FOUTNAME|=""

$!PROMPTFORTEXTSTRING |VNUM|
  INSTRUCTIONS = "Enter the number of the variable to sum:"
$!PROMPTFORTEXTSTRING |FOUTNAME|
  INSTRUCTIONS = "Enter the name of the output file containing the sum:"

$!EXTENDEDCOMMAND COMMANDPROCESSORID='extendmcr'
COMMAND='QUERY.ACTIVEZONES ZNUMS'
$!EXTENDEDCOMMAND COMMANDPROCESSORID='extendmcr'
COMMAND='QUERY.VARNAMEBYNUM |VNUM| VNAME'

$!WRITEDATASET  "|MFBD||FOUTNAME|"
 INCLUDETEXT = NO
 INCLUDEGEOM = NO
 INCLUDECUSTOMLABELS = NO
 ASSOCIATELAYOUTWITHDATAFILE = NO
 ZONELIST =  [|ZNUMS|]
 VARPOSITIONLIST =  [|VNUM|]
 BINARY = NO
 USEPOINTFORMAT = NO
 PRECISION = 14
$!RemoveVar |MFBD|
$!SYSTEM "echo Active zone |ZNUMS| > |FOUTNAME|.sum"
$!SYSTEM "echo Variable summed |VNAME|, v|VNUM|  >> |FOUTNAME|.sum"
$!SYSTEM "sort -u |FOUTNAME| | awk '{ for (i=1; i <= NF; i++) {print $i}  }' | awk '{ SUM += $1} END { print SUM }' >> |FOUTNAME|.sum"
$!SYSTEM "mv -f |FOUTNAME|.sum |FOUTNAME|"
$!SYSTEM "echo '|DATASETFNAME|' > sumvar-batch.par"
$!SYSTEM "echo |ZNUMS| >> sumvar-batch.par"
$!SYSTEM "echo |VNUM| >> sumvar-batch.par"
$!SYSTEM "echo |FOUTNAME| >> sumvar-batch.par"
$!ATTACHTEXT
  XYPOS
    {
    X = 1
    Y = 1
    }
  TEXT = "Sum of |VNAME|, v|VNUM| over active zones [|ZNUMS|] is saved on file |FOUTNAME|"
