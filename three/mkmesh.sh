#!/bin/sh

genbox << EOF
section2d.box
EOF
#
reatore2 << EOF
box
taylor
EOF
#
rm -f box.rea
rm -f taylor.rea
genmap << EOF
taylor
0.00001
EOF

