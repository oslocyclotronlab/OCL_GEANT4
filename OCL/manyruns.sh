for macrofile in sim*.mac;
do ./OCL $macrofile;
cp ../data/Edep.root ../data/${macrofile}.root;
done
