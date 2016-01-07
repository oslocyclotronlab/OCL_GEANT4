for filename in 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18;
do
root -b -q "PrintSimMacs.cpp(${filename})";
cp macfile.txt sim${filename}.mac;
done
