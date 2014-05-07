setenv PERL5LIB   /home/mikolaj/local/opt/QMPB/qmpb
setenv MEADPATH   /home/mikolaj/local
setenv QMPB       /home/mikolaj/local/opt/QMPB

#-----------------------------------------------------
set f = lysozyme

# Write QMPB input
(perl multiflex2qmpb.pl $f > $f.qmpb-in) >& multiflex2qmpb.err-out

#(perl $QMPB/tools/multiflex2qmpb.pl $f > $f.qmpb-in)>& multiflex2qmpb.err-out

# Make MEAD input
($QMPB/qmpb/qmpb.pl  $f.qmpb-in -b -inputorder > $f.qmpb-out) >& $f.err-out

# Run MEAD
cd qmpb
sh job.sh >& ../mead.err-out
cd ..

# Make GMCT input
($QMPB/qmpb/qmpb.pl  $f.qmpb-in -a -gmct -inputorder > $f.qmpb-out2)>& $f.err-out2
