#!/bin/bash
DESTDIR=$1

if [ "X$DESTDIR" == "X" ]; then
    echo "You must specify the destination directory on the command line"
    exit 1;
fi

echo -n "Copying data files..."
for file in ../data/*; do
    filename=`basename $file`
    if [ -f $file ]; then
        if [ ! -f $DESTDIR/$filename ]; then
            cp -v $file $DESTDIR
        fi
    elif [ -d $file ]; then
        if [ ! -d $DESTDIR/$filename ]; then
            cp -v -R $file $DESTDIR
        fi
    fi
done
echo "done"
