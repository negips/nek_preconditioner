#!/bin/sh

genbox << EOF
section2d.box
EOF
#
reatore2 << EOF
box
one
EOF
#
rm -f box.rea
rm -f one.rea
genmap << EOF
one
0.00001
EOF

