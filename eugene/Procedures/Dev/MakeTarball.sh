#!/bin/sh
## Sh : clean for package
echo Ce programme va supprimer des donnees dans $EUGENEDIR
echo "Voulez vous continuer ? (o/n)"
read Resp
if [ $Resp != "o" ]
then
        echo cancel ...
        exit
else
    echo on continue...
fi
cd $EUGENEDIR
make clean
./reconf

cd $EUGENEDIR/..;
find $EUGENEDIR -name ".CVS" -type d -exec rm -rf '{}' \;

cd $EUGENEDIR/..; tar czvf eugene_tarball.tar.gz eugene
echo tarball cree : $EUGENEDIR/eugene_tarball.tar.gz
