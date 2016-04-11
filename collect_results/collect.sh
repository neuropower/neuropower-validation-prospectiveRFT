BASEFOLDER="/Users/Joke/Documents/Onderzoek/Studie_4_propow/InterimPower_Results"
OUTDIR=$(echo $BASEFOLDER\tables/)

python estimation_table.py 15 80 250 'hcp' 'adaptive' 'u2' $BASEFOLDER $OUTDIR
python estimation_powpre.py 15 80 250 47 'hcp' 'adaptive' 'u2' $BASEFOLDER $OUTDIR
python estimation_powtrue.py 15 80 250 47 'hcp' 'adaptive' 'u2' $BASEFOLDER $OUTDIR
