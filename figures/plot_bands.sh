
# Get the directory the script is in, absolute path.
DIR0=$(cd $(dirname -- $0) && pwd)
FIGDIR="$DIR0/band_figures"
# Specify the common arguments.
ARGS="--ymin -1 --ymax 2.5 --zero-line --proj C,N,Cu --format png"
# Make the plots (and clean up the excess files).
echo "Plotting free bands"
sumo-bandplot $ARGS --zero-energy -2.1871 --prefix free -f "$DIR0/../HAT-Cu_free/3b.bands/vasprun.xml.gz"
rm free_band.dat
mv free_band.png "$FIGDIR"
echo "Plotting PW bands"
sumo-bandplot $ARGS --zero-energy -2.5185 --prefix PW -f "$DIR0/../HAT-Cu_PW/2b.bands/vasprun.xml.gz"
rm PW_band.dat
mv PW_band.png "$FIGDIR"
echo "Plotting Y bands"
sumo-bandplot $ARGS --zero-energy -2.4468 --prefix Y -f "$DIR0/../HAT-Cu_Y/2b.bands/vasprun.xml.gz"
rm Y_band.dat
mv Y_band.png "$FIGDIR"
rm sumo-bandplot.log
