for macrofile in sim*.mac;
do ./SingleScint $macrofile;
cp ../data/Edep.root ../data/${macrofile}.root;
done
