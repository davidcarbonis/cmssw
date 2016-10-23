#! /bin/csh
setenv RELEASE $CMSSW_VERSION

if ( ! -d /afs/cern.ch/cms/performance/tracker/activities/validation/$RELEASE/) mkdir /afs/cern.ch/cms/performance/tracker/activities/validation/$RELEASE/

setenv WWWDIR /afs/cern.ch/cms/performance/tracker/activities/validation/$RELEASE/SimHits

if ( ! -d $WWWDIR) mkdir $WWWDIR

mkdir $WWWDIR/eps
mkdir $WWWDIR/eps/TK-strips
mkdir $WWWDIR/eps/TK-strips/eloss
mkdir $WWWDIR/eps/TK-strips/Localx
mkdir $WWWDIR/eps/TK-strips/Localy
mkdir $WWWDIR/eps/TK-strips/Entryx-Exitx
mkdir $WWWDIR/eps/TK-strips/Entryy-Exity
mkdir $WWWDIR/eps/TK-strips/Entryz-Exitz
mkdir $WWWDIR/eps/TK-summary
mkdir $WWWDIR/eps/ToF

mkdir $WWWDIR/gif
mkdir $WWWDIR/gif/TK-strips
mkdir $WWWDIR/gif/TK-strips/eloss
mkdir $WWWDIR/gif/TK-strips/Localx
mkdir $WWWDIR/gif/TK-strips/Localy
mkdir $WWWDIR/gif/TK-strips/Entryx-Exitx
mkdir $WWWDIR/gif/TK-strips/Entryy-Exity
mkdir $WWWDIR/gif/TK-strips/Entryz-Exitz
mkdir $WWWDIR/gif/TK-summary
mkdir $WWWDIR/gif/ToF

echo "...Copying..."

mv plots/muon/eloss_T*_KS*.eps.gz $WWWDIR/eps/TK-strips/eloss
mv plots/muon/pos_Entryx-Exitx_T*.eps.gz $WWWDIR/eps/TK-strips/Entryx-Exitx
mv plots/muon/pos_Entryy-Exity_T*.eps.gz $WWWDIR/eps/TK-strips/Entryy-Exity
mv plots/muon/pos_Entryz-Exitz_T*.eps.gz $WWWDIR/eps/TK-strips/Entryz-Exitz
mv plots/muon/pos_Localy_T*.eps.gz $WWWDIR/eps/TK-strips/Localy
mv plots/muon/pos_Localx_T*.eps.gz $WWWDIR/eps/TK-strips/Localx

mv plots/muon/Tof.eps.gz       $WWWDIR/eps/ToF/
mv plots/muon/*summary*.eps.gz $WWWDIR/eps/TK-summary

mv plots/muon/eloss_T*_KS*.gif $WWWDIR/gif/TK-strips/eloss
mv plots/muon/pos_Entryx-Exitx_T*.gif $WWWDIR/gif/TK-strips/Entryx-Exitx
mv plots/muon/pos_Entryy-Exity_T*.gif $WWWDIR/gif/TK-strips/Entryy-Exity
mv plots/muon/pos_Entryz-Exitz_T*.gif $WWWDIR/gif/TK-strips/Entryz-Exitz
mv plots/muon/pos_Localy_T*.gif $WWWDIR/gif/TK-strips/Localy
mv plots/muon/pos_Localx_T*.gif $WWWDIR/gif/TK-strips/Localx

mv plots/muon/Tof.gif       $WWWDIR/gif/ToF/
mv plots/muon/*summary*.gif $WWWDIR/gif/TK-summary

echo "...Done..."
