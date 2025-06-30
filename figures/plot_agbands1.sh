
if [ "$1" = "-d" ]; then
    DOWNSAMPLED=1
    FILENAME="vasprun_ds.xml.gz"
    SUFFIX="_ds"
else
    DOWNSAMPLED=0
    FILENAME="vasprun.xml.gz"
    SUFFIX=""
fi
# Get the directory the script is in, absolute path.
DIR0=$(cd $(dirname -- $0) && pwd)
FIGDIR="$DIR0/band_figures"
# Specify the common arguments.
ARGS=(--ymin -1 --ymax 2.5 --zero-line --proj C,N,Cu --format png)
DOSARGS=(--config "$DIR0/sumo_colours.conf" --no-total --elements C,Cu,N --scale 3)
# Make the plots (and clean up the excess files).
echo "Plotting Ag3$SUFFIX bands"
sumo-bandplot ${ARGS[@]} ${DOSARGS[@]} --zero-energy -0.8827 --prefix "Ag3$SUFFIX" -f "$DIR0/../HAT-Cu-AgSlab/3c.bands/$FILENAME" --dos "$DIR0/../HAT-Cu-AgSlab/3b.dos/$FILENAME"
rm "Ag3${SUFFIX}_band.dat"
mv "Ag3${SUFFIX}_band.png" "$FIGDIR"
echo "Plotting Ag5$SUFFIX bands"
sumo-bandplot ${ARGS[@]} ${DOSARGS[@]} --zero-energy 0.2056 --prefix "Ag5$SUFFIX" -f "$DIR0/../HAT-Cu-AgSlab5/2c.bands/$FILENAME" --dos "$DIR0/../HAT-Cu-AgSlab5/2b.dos/$FILENAME"
rm "Ag5${SUFFIX}_band.dat"
mv "Ag5${SUFFIX}_band.png" "$FIGDIR"
rm sumo-bandplot.log
