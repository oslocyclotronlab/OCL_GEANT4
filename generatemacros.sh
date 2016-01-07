for filename in 0 18;
do
root -b -q "PrintSimMacs.cpp(${filename})";
cp macfile.txt sim${filename}.mac;
done
