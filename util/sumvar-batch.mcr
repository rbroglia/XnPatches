#!MC 1120
$!VarSet |MFBD| = './'
$!VARSET |ZNUMS|="SEDZNUMS"
$!VARSET |VNUM|="SEDVNUM"
$!VARSET |FOUTNAME|="SEDFOUTNAME"

$!READDATASET  '"|MFBD|/patch.01.p000.01.plt" "|MFBD|/patch.01.p001.01.plt" "|MFBD|/patch.01.p002.01.plt" "|MFBD|/patch.01.p003.01.plt" "|MFBD|/patch.01.p004.01.plt" "|MFBD|/patch.01.p005.01.plt" "|MFBD|/patch.01.p006.01.plt" "|MFBD|/patch.01.p007.01.plt" '
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
