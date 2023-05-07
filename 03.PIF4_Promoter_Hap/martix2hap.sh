#!/bin/bash
csvtk transpose -t coverage.martix \
| awk '{if( $70+$71+$72+$73+$74 <=3  && $49+$50+$51+$52+$53 <=3 ) print "InDel\t"$0;else if( $49+$50+$51+$52+$53 <=3 && $70+$71+$72+$73+$74 >=4 )print "DEL\t"$0;else print "Ref\t"$0 }'  | cut -f1,2 > PIF_hap.txt