cat <<EOF > 01_precedence.txt
6
0   2 3 4
1   3 4 5
EOF

cat <<EOF > 01_values.txt
7
3
-2
-2
-2
-4
EOF

../build/bin/mineflow --explicit 01_precedence.txt 01_values.txt 01_output.txt

# rm -f 01_precedence.txt 01_values.txt 01_output.txt
