for macrofile in sim*.mac;
do ./Scint $macrofile;
cp ../data/Edep.root ../data/${macrofile}.root;
done
